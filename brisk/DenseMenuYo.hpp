#include <iostream>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <omp.h>
#include <vector>
#include <cstdint>
#include <sys/resource.h>

#include "Kmers.hpp"
#include "hashing.hpp"
#include "Decycling.h"
#include "parameters.hpp"
#include "buckets.hpp"
#include "SuperKmerLight.hpp"

#include "ankerl/unordered_dense.h"

#ifndef DENSEMENUYO_H
#define DENSEMENUYO_H

#define GROGRO_THRESHOLD ((uint64_t)1 << 30)
#define bucket_overload ((uint64_t)1 << 1)



// uint64_t hash_64(uint64_t key, uint64_t mask)
// Replaced by bfc_hash_64 from hashing.hpp



template <class DATA>
class DenseMenuYo{
public:

	// Normal Buckets
	uint32_t * bucket_indexes;
	vector<Bucket<DATA> > * bucketMatrix;
	// /!\ The number of buckets is the number of distinct minimizers / 4
	uint64_t bucket_number;
	uint64_t total_number_kmers;
	uint64_t active_buckets;
	uint64_t matrix_column_number;
	// Useful variables
	kint mini_reduc_mask,m_mask;
	Parameters params;
	uint32_t enumeration_skmer_idx ;
	uint32_t enumeration_kmer_idx ;
	uint32_t largest_bucket;

	// Mutexes
	uint64_t mutex_order;
	uint64_t mutex_number;

	// Mutex
	vector<omp_lock_t> MutexData;
	vector<omp_lock_t> MutexBucket;

	// threshold kmer buckets
	ankerl::unordered_dense::map<kint, DATA> overload_kmers[bucket_overload];
	omp_lock_t lock_overload[bucket_overload];

	DATA debug_value;


	DenseMenuYo(Parameters & parameters);
	~DenseMenuYo();
	DATA * insert_kmer(kmer_full & kmer);
	DATA * get_kmer(const kmer_full & kmer);
	void insert_kmer_vector(const vector<kmer_full> & kmers, vector<bool>& newly_inserted);
	vector<DATA*> get_kmer_vector(const vector<kmer_full> & kmers) ;
	DATA * insert_kmer_no_mutex(const kmer_full & kmer, bool& newly_inserted_element,uint32_t & bucket_size, bool already_checked=false);
	DATA * get_kmer_no_mutex(const kmer_full & kmer);

	void protect_data(const kmer_full & kmer);
	void unprotect_data(const kmer_full & kmer);

	bool next(kmer_full & kmer);
	void restart_kmer_enumeration();
	void stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & largest_bucket) const;

	void print_bigest_bucket();

private:
	friend class BriskWriter;

	/** Translate a bucket into an internal map and remove the bucket
	 */
	void bucket_to_map(uint64_t small_minimizer);

	// Buffers
	kint buffered_minimizer;
	kmer_full * buffered_kmer;
	
	bool enumeration_started;
	typename ankerl::unordered_dense::map<kint, DATA>::iterator overload_iter;
	uint32_t current_overload;
	uint64_t current_minimizer;
};



template <class DATA>
DenseMenuYo<DATA>::DenseMenuYo(Parameters & parameters): params( parameters ) {
	// If m - reduc > 16, not enougth space for bucket idxs ! (because of uint32_t)
	this->mini_reduc_mask = ((kint)1 << (2 * params.b)) - 1;
	this->m_mask = ((kint)1 << (2 * params.m)) - 1;

	// Create the mutexes
	mutex_order = min((uint8_t)10,params.b);
	mutex_number = ((uint64_t)1)<<(2*mutex_order); /* DO NOT MODIFY, to increase/decrease modify mutex_order*/
	// Init the mutexes
	MutexData.resize(mutex_number);
	MutexBucket.resize(mutex_number);
	for (uint64_t i(0); i < mutex_number; ++i) {
		omp_init_lock(&MutexData[i]);
		omp_init_lock(&MutexBucket[i]);
	}
	total_number_kmers=0;
	active_buckets=0;
	largest_bucket=0;

	bucket_number=(uint64_t)1<<(2*params.b);
	matrix_column_number = bucket_number / mutex_number;
	bucket_indexes = new uint32_t[bucket_number];
	memset(bucket_indexes, 0, sizeof(uint32_t)*bucket_number);
	bucketMatrix = new vector<Bucket<DATA> >[mutex_number];
	// for(uint32_t i(0);i<mutex_number;++i){
	// 	bucketMatrix[i].reserve(matrix_column_number);
	// }

	enumeration_started = false;
	current_minimizer = 0;
	current_overload=0;
	enumeration_skmer_idx = 0;
	enumeration_kmer_idx = 0;
	overload_iter=overload_kmers[current_overload].begin();
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
void DenseMenuYo<DATA>::bucket_to_map(uint64_t small_minimizer) {
	// Get the bucket
	uint32_t mutex_idx = get_mutex(small_minimizer);
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);
	uint32_t idx = bucket_indexes[matrix_idx];
	uint64_t kmer_id(0);

	Bucket<DATA> &bucket = bucketMatrix[mutex_idx][idx-1];

	// Init values for copy
	kmer_full kmer;
	DATA * data;

	// Enumerate and save all values
	omp_set_lock(&lock_overload[small_minimizer%bucket_overload]);
	// bucket.init_enumerator();
	uint32_t enumeration_skmer_idx = 0;
	uint32_t enumeration_kmer_idx = 0;
	while (bucket.has_next_kmer(enumeration_skmer_idx,enumeration_kmer_idx)) {
		bucket.next_kmer(kmer, small_minimizer,enumeration_skmer_idx,enumeration_kmer_idx);
		kmer.minimizer_idx -= (params.m_reduc + 1)/2;
		kmer.compute_mini(params.m);
		data = bucket.data_reserved_memory+kmer_id*sizeof(DATA);
		kmer.unhash_kmer_minimizer(params.m);
		kint local_kmer = kmer.kmer_s;
		this->overload_kmers[small_minimizer%bucket_overload][local_kmer]=*data;
		kmer_id++;
	}
	omp_unset_lock(&lock_overload[small_minimizer%bucket_overload]);

	// Free some memory from the bucket
	bucket.clear();
}



//NON c'est NON
template <class DATA>
DATA * DenseMenuYo<DATA>::insert_kmer(kmer_full & kmer) {
	DATA * prev_val = this->get_kmer(kmer);
	if (prev_val != NULL) {
		return prev_val;
	}

	// Transform the super minimizer to the used minimizer
	// Remove the minimizer suffix
	uint64_t small_minimizer = kmer.minimizer >> (2 * ((this->params.m_reduc + 1) / 2));
	// Remove the minimizer prefix
	small_minimizer &= this->mini_reduc_mask;

	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Works because m - reduc <= 16
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	// omp_set_lock(&MutexBucket[mutex_idx]);
	uint32_t idx = bucket_indexes[matrix_idx];

	// GROGRO bucket check
	if (bucketMatrix[mutex_idx][idx-1].sorted_size<0) {
		omp_set_lock(&lock_overload[small_minimizer%bucket_overload]);
		kint local_kmer = kmer.get_unhash_kmer_value(params.m);
		overload_kmers[small_minimizer%bucket_overload][local_kmer] = DATA();
		DATA* el=&overload_kmers[small_minimizer%bucket_overload][local_kmer];
		// bucketMatrix[mutex_idx][idx-1].nb_kmers += 1;
		omp_unset_lock(&lock_overload[small_minimizer%bucket_overload]);
		return el;
	}

	// Create the bucket if not already existing
	if (idx == 0) {
		bucketMatrix[mutex_idx].emplace_back(&params);
		bucket_indexes[matrix_idx] = bucketMatrix[mutex_idx].size();
		idx = bucket_indexes[matrix_idx];
		active_buckets++;
	}
	// Insert the kmer in the right bucket
	DATA * value = bucketMatrix[mutex_idx][idx-1].insert_kmer(kmer);

	// Bucket now too big to hold
	if ((not (bucketMatrix[mutex_idx][idx-1].sorted_size<0)) and bucketMatrix[mutex_idx][idx-1].skml.size() >= GROGRO_THRESHOLD) {
		this->bucket_to_map(small_minimizer);
	}

	return value;
}



template <class DATA>
void DenseMenuYo<DATA>::insert_kmer_vector(const vector<kmer_full> & kmers,vector<bool>& newly_inserted) {
	vector<DATA *> result = this->get_kmer_vector(kmers);
	uint32_t bucket_size = 0;
	for(uint i(0);i<kmers.size(); ++i){
		if(result[i]==NULL){
		// if(true){
			bool new_insert = false;
			insert_kmer_no_mutex(kmers[i],new_insert,bucket_size,true);
			newly_inserted.push_back(new_insert);
		}else{
			newly_inserted.push_back(false);
		}
	}

	if (bucket_size > largest_bucket) {
		largest_bucket = bucket_size;
	}

	return;
}



template <class DATA>
DATA * DenseMenuYo<DATA>::insert_kmer_no_mutex(const kmer_full & kmer,bool& newly_inserted,uint32_t & bucket_size, bool already_checked) {
	if(not already_checked){
		DATA * prev_val = this->get_kmer_no_mutex(kmer);
		if (prev_val != NULL) {
			newly_inserted=false;
			return prev_val;
		}
	}
	newly_inserted=true;
	// Transform the super minimizer to the used minimizer

	// Remove the minimizer suffix
	uint64_t small_minimizer = kmer.minimizer >> (2 * ((this->params.m_reduc + 1) / 2));
	// Remove the minimizer prefix
	small_minimizer &= this->mini_reduc_mask;
	
	// cout << "insert small mini " << kmer2str(small_minimizer, params.b) << endl;

	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Works because m - reduc <= 16
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);
	uint32_t idx = bucket_indexes[matrix_idx];
	
	// Create the bucket if not already existing
	if (idx == 0) {
		// cout<<"emplaceback"<<endl;
		bucketMatrix[mutex_idx].emplace_back();
		// cout<<"ENDEB"<<endl;
		bucket_indexes[matrix_idx] = bucketMatrix[mutex_idx].size();
		idx = bucket_indexes[matrix_idx];
		active_buckets++;
	}

	// Insert the kmer in the right bucket
	DATA * value = bucketMatrix[mutex_idx][idx-1].insert_kmer(kmer);
	#pragma omp atomic
	total_number_kmers++;

	bucket_size = bucketMatrix[mutex_idx][idx-1].skml.size();

	return value;
}



template <class DATA>
DATA * DenseMenuYo<DATA>::get_kmer(const kmer_full & kmer) {
	// Remove the minimizer suffix
	uint64_t small_minimizer = kmer.minimizer >> (2 * ((this->params.m_reduc + 1) / 2));
	// Remove the minimizer prefix
	small_minimizer &= this->mini_reduc_mask;


	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Normal kmer
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	// omp_set_lock(&MutexBucket[mutex_idx]);
	uint32_t idx = bucket_indexes[matrix_idx];
	// No bucket
	if (idx == 0) {
		// cout << "PAGLOP idx " << kmer2str(small_minimizer, params.b) << endl;
		return (DATA *)NULL;
	}

	// GROGRO bucket check
	if (bucketMatrix[mutex_idx][idx-1].sorted_size<0) {
		kint local_kmer = kmer.get_unhash_kmer_value(params.m);
		omp_set_lock(&lock_overload[small_minimizer%bucket_overload]);
		if (overload_kmers[small_minimizer%bucket_overload].count(local_kmer) != 0) {
			DATA* el=&(overload_kmers[small_minimizer%bucket_overload][local_kmer]);
			omp_unset_lock(&lock_overload[small_minimizer%bucket_overload]);
			return el;
		} else {
			omp_unset_lock(&lock_overload[small_minimizer%bucket_overload]);
			return (DATA *)NULL;
		}
	}

	// Looks into the bucket for the right kmer
	DATA * value = bucketMatrix[mutex_idx][idx-1].find_kmer(kmer);

	return value;
}



template <class DATA>
vector<DATA *> DenseMenuYo<DATA>::get_kmer_vector(const vector<kmer_full> & kmers)  {
	vector<DATA *> result(kmers.size(),NULL);
	//WE ASSUME HERE THAT ALL KMER HAVE THE SAME MINIMIZER
    // Remove the minimizer suffix
	uint64_t small_minimizer = kmers[0].minimizer >> (2 * ((this->params.m_reduc + 1) / 2));
	// Remove the minimizer prefix
	small_minimizer &= this->mini_reduc_mask;

	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Normal kmer
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	uint32_t idx = bucket_indexes[matrix_idx];
	// No bucket
	if (idx == 0) {
		return result;
	}

	// On GROGRO bucket
	if (bucketMatrix[mutex_idx][idx-1].sorted_size<0) {
		omp_set_lock(&lock_overload[small_minimizer%bucket_overload]);
		for(uint i(0);i<kmers.size();++i){
			kint local_kmer=kmers[i].get_unhash_kmer_value(params.m);
			if (overload_kmers[small_minimizer%bucket_overload].count(local_kmer) != 0) {
				result[i] =&overload_kmers[small_minimizer%bucket_overload][local_kmer];
			}else{
			}
			
		}
		omp_unset_lock(&lock_overload[small_minimizer%bucket_overload]);
		return result;
	}
	// Looks into the bucket for the right kmer
	// cout << "small mini " << kmer2str(small_minimizer, params.b) << endl;
	result = bucketMatrix[mutex_idx][idx-1].find_kmer_vector(kmers);
	return result;
}



template <class DATA>
DATA * DenseMenuYo<DATA>::get_kmer_no_mutex(const kmer_full & kmer) {	
    // Remove the minimizer suffix
	uint64_t small_minimizer = kmer.minimizer >> (2 * ((this->params.m_reduc + 1) / 2));
	// Remove the minimizer prefix
	small_minimizer &= this->mini_reduc_mask;

	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Normal kmer
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	uint32_t idx = bucket_indexes[matrix_idx];
	// No bucket
	if (idx == 0) {
		//cout << "PAGLOP idx " << kmer2str(small_minimizer, params.b) << endl; cin.get();
		return (DATA *)NULL;
	}

	DATA * value;
	// Looks into the bucket for the right kmer
	value = bucketMatrix[mutex_idx][idx-1].find_kmer(kmer);
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
    // Remove the minimizer suffix
	uint64_t small_minimizer = kmer.minimizer >> (2 * ((this->params.m_reduc + 1) / 2));
	// Remove the minimizer prefix
	small_minimizer &= this->mini_reduc_mask;

	uint32_t mutex_idx = get_mutex(small_minimizer);
	omp_set_lock(&MutexData[mutex_idx]);
}



template<class DATA>
void DenseMenuYo<DATA>::unprotect_data(const kmer_full & kmer) {
    // Remove the minimizer suffix
	uint64_t small_minimizer = kmer.minimizer >> (2 * ((this->params.m_reduc + 1) / 2));
	// Remove the minimizer prefix
	small_minimizer &= this->mini_reduc_mask;

	uint32_t mutex_idx = get_mutex(small_minimizer);
	omp_unset_lock(&MutexData[mutex_idx]);
}



template<class DATA>
bool DenseMenuYo<DATA>::next(kmer_full & kmer) {
	if (not enumeration_started) {
		current_overload=0;
		overload_iter=overload_kmers[current_overload].begin();
		enumeration_started = true;
	}
	// cout<<"next0"<<endl;

	uint32_t mutex_idx = get_mutex(current_minimizer);
	uint32_t column_idx = get_column(current_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);
	uint32_t idx = bucket_indexes[matrix_idx];
	// Select minimizer
	while (current_minimizer<bucket_number and idx == 0) {
		mutex_idx = get_mutex(current_minimizer);
		column_idx = get_column(current_minimizer);
		matrix_idx = get_matrix_position(mutex_idx, column_idx);
		idx = bucket_indexes[matrix_idx];
		// cout<<"next1"<<endl;

		// Already enumerated kmers
		if (idx != 0 and bucketMatrix[mutex_idx][idx-1].sorted_size<0) // Maybe a bug here
			idx = 0;

		if (idx == 0){
			current_minimizer += 1;
			enumeration_skmer_idx = 0;
			enumeration_kmer_idx = 0;
		}
	}

	// Enumeration is over
	if (current_minimizer >= bucket_number){
		// cout<<"next2"<<endl;
		return false;
	}

	if (not bucketMatrix[mutex_idx][idx-1].has_next_kmer(enumeration_skmer_idx,enumeration_kmer_idx)) {
		// cout<<"next3"<<endl;
		current_minimizer += 1;
		enumeration_skmer_idx = 0;
		enumeration_kmer_idx = 0;
		return this->next(kmer);
	}
	// cout<<"next4"<<endl;
	bucketMatrix[mutex_idx][idx-1].next_kmer(kmer, current_minimizer,enumeration_skmer_idx,enumeration_kmer_idx);
	kmer.minimizer_idx -= (params.m_reduc + 1)/2;
	kmer.compute_mini(params.m);

	return true;
}



template<class DATA>
void DenseMenuYo<DATA>::restart_kmer_enumeration() {
	current_minimizer = 0;
	enumeration_skmer_idx = 0;
	enumeration_kmer_idx = 0;
}


static inline uint32_t mylog2(const uint32_t x) {
  uint32_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}



template<class DATA>
void DenseMenuYo<DATA>::stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & largest_bucket) const {
	nb_buckets = 0;
	nb_kmers = 0;
	nb_skmers = 0;
	largest_bucket=0;

	for (uint64_t mini=0 ; mini<bucket_number ; mini++) {
		uint32_t mutex_idx = get_mutex(mini);
		uint32_t column_idx = get_column(mini);
		uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);
		uint32_t idx = bucket_indexes[matrix_idx];

		if (idx != 0) {
			#pragma omp atomic
			nb_buckets += 1;
			#pragma omp atomic
			nb_skmers += bucketMatrix[mutex_idx][idx-1].skml.size();
			if(largest_bucket<bucketMatrix[mutex_idx][idx-1].skml.size()){
				#pragma omp critical (max)
				largest_bucket=bucketMatrix[mutex_idx][idx-1].skml.size();
			}
		}
	}
	nb_kmers = total_number_kmers;
}



#endif
