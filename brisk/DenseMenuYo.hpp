#include <iostream>
#include <cstdint>
#include <cstring>
#include <map>
#include <omp.h>
#include <vector>

#include <sys/resource.h>

#include "Kmers.hpp"
#include "buckets.hpp"
#include "SuperKmerLight.hpp"
#include "parameters.hpp"

#ifndef DENSEMENUYO_H
#define DENSEMENUYO_H



template <class DATA>
class DenseMenuYo{
public:

	// Normal Buckets
	uint32_t * bucket_indexes;
	vector<Bucket<DATA> > * bucketMatrix;
	// /!\ The number of bucket is the number of distinct minimizer / 4
	uint64_t bucket_number;
	uint64_t matrix_column_number;
	// Usefull variables
	kint mini_reduc_mask;
	Parameters params;
	

	// Mutexes
	uint64_t mutex_order;
	uint64_t mutex_number;

	// Mutex
	vector<omp_lock_t> MutexData;
	vector<omp_lock_t> MutexBucket;
	omp_lock_t multi_lock;

	// Cursed kmer buckets
	robin_hood::unordered_map<kint, DATA> cursed_kmers;

	DATA debug_value;


	DenseMenuYo(Parameters & parameters);
	~DenseMenuYo();
	DATA * insert_kmer(kmer_full & kmer);
	DATA * get_kmer(kmer_full & kmer);
	vector<DATA*> insert_kmer_vector(vector<kmer_full> & kmers, vector<bool>& newly_inserted);
	vector<DATA*> get_kmer_vector(vector<kmer_full> & kmers);
	DATA * insert_kmer_no_mutex(kmer_full & kmer, bool& newly_inserted_element,bool already_checked=false);
	DATA * get_kmer_no_mutex(kmer_full & kmer);

	void protect_data(const kmer_full & kmer);
	void unprotect_data(const kmer_full & kmer);

	bool next(kmer_full & kmer);
	void restart_kmer_enumeration();
	void stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & nb_cursed, uint64_t & largest_bucket) const;

	void print_bigest_bucket();

private:
	friend class BriskWriter;
	
	

	

	// Buffers
	kint buffered_minimizer;
	kmer_full * buffered_kmer;
	
	// uint64_t getMemorySelfMaxUsed();
	bool enumeration_started;
	typename robin_hood::unordered_map<kint, DATA>::iterator cursed_iter;
	uint32_t current_minimizer;
};


template <class DATA>
DenseMenuYo<DATA>::DenseMenuYo(Parameters & parameters): params( parameters ) {
	// If m - reduc > 16, not enougth space for bucket idxs ! (because of uint32_t)
	assert(params.m_small <= 16);
	mini_reduc_mask = ((kint)1 << (2 * params.m_small)) - 1;

	// Create the mutexes
	mutex_order = min((uint8_t)6,params.m_small);
	mutex_number = ((uint64_t)1)<<(2*mutex_order); /* DO NOT MODIFY, to increase/decrease modify mutex_order*/
	// Init the mutexes
	MutexData.resize(mutex_number);
	MutexBucket.resize(mutex_number);
	for (uint64_t i(0); i < mutex_number; ++i) {
		omp_init_lock(&MutexData[i]);
		omp_init_lock(&MutexBucket[i]);
	}
	omp_init_lock(&multi_lock);

	bucket_number=1<<(2*params.m_small);
	matrix_column_number = bucket_number / mutex_number;
	bucket_indexes = new uint32_t[bucket_number];
	memset(bucket_indexes, 0, sizeof(uint32_t)*bucket_number);
	bucketMatrix = new vector<Bucket<DATA> >[mutex_number];

	enumeration_started = false;
	cursed_iter = cursed_kmers.begin();
	current_minimizer = 0;
}


template <class DATA>
DenseMenuYo<DATA>::~DenseMenuYo() {
	delete[] bucket_indexes;
	delete[] bucketMatrix;
}


#define get_mutex(mini) (mini%mutex_number)
#define get_column(mini) (mini>>(2*mutex_order))
#define get_matrix_position(row_idx, col_idx) (row_idx * matrix_column_number + col_idx)


template <class DATA>
DATA * DenseMenuYo<DATA>::insert_kmer(kmer_full & kmer) {
	DATA * prev_val = this->get_kmer(kmer);
	if (prev_val != NULL) {
		return prev_val;
	}

	// Cursed kmers
	if (kmer.multi_mini) {
		omp_set_lock(&multi_lock);
		cursed_kmers[kmer.kmer_s] = DATA();
		omp_unset_lock(&multi_lock);
		// return &debug_value;
		return &(cursed_kmers[kmer.kmer_s]);
	}

	// Transform the super minimizer to the used minimizer
	uint64_t small_minimizer = (uint32_t)(kmer.minimizer & mini_reduc_mask);
	// kmer.minimizer_idx += params.m_reduc;
	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Works because m - reduc <= 16
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	omp_set_lock(&MutexBucket[mutex_idx]);
	uint32_t idx = bucket_indexes[matrix_idx];
	// Create the bucket if not already existing
	if (idx == 0) {
		bucketMatrix[mutex_idx].emplace_back(&params);
		bucket_indexes[matrix_idx] = bucketMatrix[mutex_idx].size();
		idx = bucket_indexes[matrix_idx];
	}
	// Insert the kmer in the right bucket
	DATA * value = bucketMatrix[mutex_idx][idx-1].insert_kmer(kmer);
	omp_unset_lock(&MutexBucket[mutex_idx]);

	// kmer.minimizer_idx -= params.m_reduc;
	return value;
}


template <class DATA>
vector<DATA *> DenseMenuYo<DATA>::insert_kmer_vector(vector<kmer_full> & kmers,vector<bool>& newly_inserted) {
	vector<DATA *> result = this->get_kmer_vector(kmers);
	for(uint i(0);i<kmers.size(); ++i){
		if(result[i]==NULL){
			newly_inserted.push_back(true);
			bool dog;
			insert_kmer_no_mutex(kmers[i],dog,true);
		}else{
			newly_inserted.push_back(false);
		}
	}
	// uint64_t small_minimizer = (uint32_t)(kmers[0b11].minimizer & mini_reduc_mask);
	// // kmer.minimizer_idx += params.m_reduc;
	// uint32_t mutex_idx = get_mutex(small_minimizer);

	// // Works because m - reduc <= 16
	// uint32_t column_idx = get_column(small_minimizer);
	// uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	// uint32_t idx = bucket_indexes[matrix_idx];
	// // Create the bucket if not already existing
	// if (idx == 0) {
	// 	bucketMatrix[mutex_idx].emplace_back(&params);
	// 	bucket_indexes[matrix_idx] = bucketMatrix[mutex_idx].size();
	// 	idx = bucket_indexes[matrix_idx];
	// }
	// for(uint i(0);i<kmers.size(); ++i){
	// 	if(result[i]==NULL){
	// 		newly_inserted.push_back(true);
	// 		if (kmers[i].multi_mini) {
	// 			omp_set_lock(&multi_lock);
	// 			cursed_kmers[kmers[i].kmer_s] = DATA();
	// 			omp_unset_lock(&multi_lock);
	// 			result[i]=&(cursed_kmers[kmers[i].kmer_s]);
	// 		}else{
	// 			result[i] = bucketMatrix[mutex_idx][idx-1].insert_kmer(kmers[i]);
	// 		}
	// 	}else{
	// 		newly_inserted.push_back(false);
	// 	}
	// }
	return result;
}


template <class DATA>
DATA * DenseMenuYo<DATA>::insert_kmer_no_mutex(kmer_full & kmer,bool& newly_inserted, bool already_checked) {
	if(not already_checked){
		DATA * prev_val = this->get_kmer_no_mutex(kmer);
		if (prev_val != NULL) {
			newly_inserted=false;
			return prev_val;
		}
	}
	newly_inserted=true;
	// Cursed kmers
	if (kmer.multi_mini) {
		omp_set_lock(&multi_lock);
		cursed_kmers[kmer.kmer_s] = DATA();
		omp_unset_lock(&multi_lock);
		// return &debug_value;
		return &(cursed_kmers[kmer.kmer_s]);
	}
	// Transform the super minimizer to the used minimizer
	uint64_t small_minimizer = (uint32_t)(kmer.minimizer & mini_reduc_mask);
	// kmer.minimizer_idx += params.m_reduc;
	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Works because m - reduc <= 16
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	
	uint32_t idx = bucket_indexes[matrix_idx];
	// Create the bucket if not already existing
	if (idx == 0) {
		bucketMatrix[mutex_idx].emplace_back(&params);
		bucket_indexes[matrix_idx] = bucketMatrix[mutex_idx].size();
		idx = bucket_indexes[matrix_idx];
	}
	// Insert the kmer in the right bucket
	DATA * value = bucketMatrix[mutex_idx][idx-1].insert_kmer(kmer);

	return value;
}



template <class DATA>
DATA * DenseMenuYo<DATA>::get_kmer(kmer_full & kmer) {
	// Cursed kmers
	if (kmer.multi_mini) {
		omp_set_lock(&multi_lock);
		if (cursed_kmers.count(kmer.kmer_s) != 0) {
			omp_unset_lock(&multi_lock);
			return &(cursed_kmers[kmer.kmer_s]);
		} else {
			omp_unset_lock(&multi_lock);
			return (DATA *)NULL;
		}
	}
	
	uint64_t small_minimizer = (uint32_t)(kmer.minimizer & mini_reduc_mask);
	// kmer.minimizer_idx += params.m_reduc;
	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Normal kmer
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	omp_set_lock(&MutexBucket[mutex_idx]);
	uint32_t idx = bucket_indexes[matrix_idx];
	// No bucket
	if (idx == 0) {
		omp_unset_lock(&MutexBucket[mutex_idx]);
		// kmer.minimizer_idx -= params.m_reduc;
		return (DATA *)NULL;
	}

	DATA * value;
	// Looks into the bucket for the right kmer
	value = bucketMatrix[mutex_idx][idx-1].find_kmer(kmer);

	omp_unset_lock(&MutexBucket[mutex_idx]);	
	// kmer.minimizer_idx -= params.m_reduc;
	return value;
}


template <class DATA>
vector<DATA *> DenseMenuYo<DATA>::get_kmer_vector(vector<kmer_full> & kmers) {
	// cout<<"get_kmer_vector"<<endl;
	vector<DATA *> result(kmers.size(),NULL);
	//WE ASSUME HERE THAT ALL KMERS are equally cursed
	if (kmers[0].multi_mini) {
		// cout<<"cvursedf"<<endl;
		omp_set_lock(&multi_lock);
		for(uint i(0);i<kmers.size();++i){
			// if(kmers[i].multi_mini==false){
			// 	cout<<"WTF MAN"<<endl;
			// 	cin.get();
			// }
			if (cursed_kmers.count(kmers[i].kmer_s) != 0) {
				result[i] =&cursed_kmers[kmers[i].kmer_s];
			}
		}
		omp_unset_lock(&multi_lock);
		return result;
	}
	
	//WE ASSUME HERE THAT ALL KMER HAVE THE SAME MINIMIZER
	uint64_t small_minimizer = (uint32_t)(kmers[0].minimizer & mini_reduc_mask);
	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Normal kmer
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	uint32_t idx = bucket_indexes[matrix_idx];
	// No bucket
	if (idx == 0) {
		// cout<<"OUT get_kmer_vector OUT "<<endl;

		return result;
	}

	// Looks into the bucket for the right kmer
	result = bucketMatrix[mutex_idx][idx-1].find_kmer_vector(kmers);
	// for(uint i(0);i<result.size();i++) {
	// 	if((result[i])!=bucketMatrix[mutex_idx][idx-1].find_kmer(kmers[i])){
	// 		cout<<"FIND WTF"<<endl;
	// 		if((result[i])==NULL){cout<<"RNULL"<<endl;}
	// 		if(bucketMatrix[mutex_idx][idx-1].find_kmer(kmers[i])==NULL){cout<<"F  NULL"<<endl;}
	// 		// cin.get();
	// 	}
	// 	if((result[i])!=get_kmer_no_mutex(kmers[i])){
			
	// 		cout<<"FIND WTF2"<<endl;
	// 		cin.get();
	// 	}else{
	// 		// cout<<i<<" ok"<<endl;
	// 	}
	// }
	// cout<<"OUT get_kmer_vector OUT "<<endl;
	return result;
}


template <class DATA>
DATA * DenseMenuYo<DATA>::get_kmer_no_mutex(kmer_full & kmer) {
	
	// Cursed kmers
	if (kmer.multi_mini) {
		// cout<<"cursed"<<endl;
		omp_set_lock(&multi_lock);
		if (cursed_kmers.count(kmer.kmer_s) != 0) {
			omp_unset_lock(&multi_lock);
			return &(cursed_kmers[kmer.kmer_s]);
		} else {
			omp_unset_lock(&multi_lock);
			return (DATA *)NULL;
		}
	}
	
	uint64_t small_minimizer = (uint32_t)(kmer.minimizer & mini_reduc_mask);
	// kmer.minimizer_idx += params.m_reduc;
	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Normal kmer
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	uint32_t idx = bucket_indexes[matrix_idx];
	// No bucket
	if (idx == 0) {
		// kmer.minimizer_idx -= params.m_reduc;
		return (DATA *)NULL;
	}

	DATA * value;
	// Looks into the bucket for the right kmer
	value = bucketMatrix[mutex_idx][idx-1].find_kmer(kmer);
	// kmer.minimizer_idx -= params.m_reduc;
	return value;
}


template<class DATA>
void DenseMenuYo<DATA>::print_bigest_bucket() {
	uint max_size = 0;
	uint max_mutex = 0;
	uint max_idx = 0;
	for (uint mutex=0 ; mutex<mutex_number ; mutex++) {
		for (uint i=0 ; i<bucketMatrix[mutex].size() ; i++) {
			if (max_size < bucketMatrix[mutex][i].skml.size()) {
				max_size = bucketMatrix[mutex][i].skml.size();
				max_mutex = mutex;
				max_idx = i;
			}
		}
	}

	bucketMatrix[max_mutex][max_idx].print();
}


template<class DATA>
void DenseMenuYo<DATA>::protect_data(const kmer_full & kmer) {
	uint64_t small_minimizer = (uint32_t)(kmer.minimizer & mini_reduc_mask);
	uint32_t mutex_idx = get_mutex(small_minimizer);

	omp_set_lock(&MutexData[mutex_idx]);
}


template<class DATA>
void DenseMenuYo<DATA>::unprotect_data(const kmer_full & kmer) {
	uint64_t small_minimizer = (uint32_t)(kmer.minimizer & mini_reduc_mask);
	uint32_t mutex_idx = get_mutex(small_minimizer);

	omp_unset_lock(&MutexData[mutex_idx]);
}


template<class DATA>
bool DenseMenuYo<DATA>::next(kmer_full & kmer) {
	if (not enumeration_started) {
		cursed_iter = cursed_kmers.begin();
		enumeration_started = true;
	}
	// Cursed kmers
	if (cursed_iter != cursed_kmers.end()) {
		kmer.kmer_s = cursed_iter->first;
		kmer.multi_mini = true;

		cursed_iter = std::next(cursed_iter);
		return true;
	}
	
	uint32_t mutex_idx = get_mutex(current_minimizer);
	uint32_t column_idx = get_column(current_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);
	uint32_t idx = bucket_indexes[matrix_idx];
	// Select minimizer
	while (current_minimizer<bucket_number and idx == 0) {
		current_minimizer += 1;
		mutex_idx = get_mutex(current_minimizer);
		column_idx = get_column(current_minimizer);
		matrix_idx = get_matrix_position(mutex_idx, column_idx);
		idx = bucket_indexes[matrix_idx];
	}

	// Enumeration is over
	if (current_minimizer >= bucket_number)
		return false;

	if (not bucketMatrix[mutex_idx][idx-1].has_next_kmer()) {
		current_minimizer += 1;
		return this->next(kmer);
	}
	
	bucketMatrix[mutex_idx][idx-1].next_kmer(kmer, current_minimizer);
	// kmer.minimizer_idx -= params.m_reduc;
	kmer.compute_mini(params.m);

	return true;
}


template<class DATA>
void DenseMenuYo<DATA>::restart_kmer_enumeration() {
	current_minimizer = 0;
}

template<class DATA>
void DenseMenuYo<DATA>::stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & nb_cursed, uint64_t & largest_bucket) const {
	nb_buckets = 0;
	nb_kmers = 0;
	nb_skmers = 0;
	largest_bucket=0;

	//#pragma omp parallel for
	for (uint32_t mini=0 ; mini<bucket_number ; mini++) {
		uint32_t mutex_idx = get_mutex(mini);
		uint32_t column_idx = get_column(mini);
		uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);
		uint32_t idx = bucket_indexes[matrix_idx];

		if (idx != 0) {
			#pragma omp atomic
			nb_buckets += 1;
			#pragma omp atomic
			nb_skmers += bucketMatrix[mutex_idx][idx-1].skml.size();
			#pragma omp atomic
			nb_kmers += bucketMatrix[mutex_idx][idx-1].nb_kmers;
			if(largest_bucket<bucketMatrix[mutex_idx][idx-1].skml.size()){
				#pragma omp critical (max)
				largest_bucket=bucketMatrix[mutex_idx][idx-1].skml.size();
			}
		}
	}

	nb_cursed = cursed_kmers.size();
	nb_kmers += nb_cursed;
}



#endif
