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

#define GROGRO_THREASHOLD ((uint64_t)1 << 3)
#define bucket_overload ((uint64_t)1 << 12)



uint64_t hash_64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}



// The inversion of hash_64(). Modified from <https://naml.us/blog/tag/invertible>
uint64_t hash_64i(uint64_t key, uint64_t mask)
{
	uint64_t tmp;

	// Invert key = key + (key << 31)
	tmp = (key - (key << 31));
	key = (key - (tmp << 31)) & mask;

	// Invert key = key ^ (key >> 28)
	tmp = key ^ key >> 28;
	key = key ^ tmp >> 28;

	// Invert key *= 21
	key = (key * 14933078535860113213ull) & mask;

	// Invert key = key ^ (key >> 14)
	tmp = key ^ key >> 14;
	tmp = key ^ tmp >> 14;
	tmp = key ^ tmp >> 14;
	key = key ^ tmp >> 14;

	// Invert key *= 265
	key = (key * 15244667743933553977ull) & mask;

	// Invert key = key ^ (key >> 24)
	tmp = key ^ key >> 24;
	key = key ^ tmp >> 24;

	// Invert key = (~key) + (key << 21)
	tmp = ~key;
	tmp = ~(key - (tmp << 21));
	tmp = ~(key - (tmp << 21));
	key = ~(key - (tmp << 21)) & mask;

	return key;
}



template <class DATA>
class DenseMenuYo{
public:

	// Normal Buckets
	uint32_t * bucket_indexes;
	vector<Bucket<DATA> > * bucketMatrix;
	// /!\ The number of bucket is the number of distinct minimizer / 4
	uint64_t bucket_number;
	uint64_t total_number_kmers;
	uint64_t active_buckets;
	uint64_t matrix_column_number;
	// Usefull variables
	kint mini_reduc_mask,m_mask;
	Parameters params;
	

	// Mutexes
	uint64_t mutex_order;
	uint64_t mutex_number;

	// Mutex
	vector<omp_lock_t> MutexData;
	vector<omp_lock_t> MutexBucket;
	omp_lock_t multi_lock;

	// Cursed kmer buckets
	robin_hood::unordered_node_map<kint, DATA> cursed_kmers;
	robin_hood::unordered_node_map<kint, DATA> overload_kmers[bucket_overload];
	// tsl::sparse_map<kint, DATA> cursed_kmers;
	// tsl::sparse_map<kint, DATA> overload_kmers[bucket_overload];
	omp_lock_t lock_overload[bucket_overload];

	DATA debug_value;


	DenseMenuYo(Parameters & parameters);
	~DenseMenuYo();
	DATA * insert_kmer(kmer_full & kmer);
	DATA * get_kmer(kmer_full & kmer);
	void insert_kmer_vector(const vector<kmer_full> & kmers, vector<bool>& newly_inserted);
	vector<DATA*> get_kmer_vector(const vector<kmer_full> & kmers) ;
	DATA * insert_kmer_no_mutex(const kmer_full & kmer, bool& newly_inserted_element,bool already_checked=false);
	DATA * get_kmer_no_mutex(const kmer_full & kmer);

	void protect_data(const kmer_full & kmer);
	void unprotect_data(const kmer_full & kmer);

	bool next(kmer_full & kmer);
	void restart_kmer_enumeration();
	void stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & nb_cursed, uint64_t & largest_bucket) const;

	void print_bigest_bucket();

private:
	friend class BriskWriter;

	/** Translate a bucket into an internal map and remove the bucket
	 */
	void bucket_to_map(uint64_t small_minimizer);

	// Buffers
	kint buffered_minimizer;
	kmer_full * buffered_kmer;
	
	// uint64_t getMemorySelfMaxUsed();
	bool enumeration_started;
	typename robin_hood::unordered_node_map<kint, DATA>::iterator cursed_iter;
	typename robin_hood::unordered_node_map<kint, DATA>::iterator overload_iter;
	// typename tsl::sparse_map<kint, DATA>::iterator cursed_iter;
	// typename tsl::sparse_map<kint, DATA>::iterator overload_iter;
	uint32_t current_overload;
	uint64_t current_minimizer;
};



template <class DATA>
DenseMenuYo<DATA>::DenseMenuYo(Parameters & parameters): params( parameters ) {
	// If m - reduc > 16, not enougth space for bucket idxs ! (because of uint32_t)
	mini_reduc_mask = ((kint)1 << (2 * params.m_small)) - 1;
	mini_reduc_mask<<=((params.m-params.m_small));
	m_mask = ((kint)1 << (2 * params.m)) - 1;

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
	for (uint64_t i(0); i < bucket_overload; ++i) {
		omp_init_lock(&lock_overload[i]);
	}
	total_number_kmers=0;
	active_buckets=0;

	bucket_number=(uint64_t)1<<(2*params.m_small);
	matrix_column_number = bucket_number / mutex_number;
	bucket_indexes = new uint32_t[bucket_number];
	memset(bucket_indexes, 0, sizeof(uint32_t)*bucket_number);
	bucketMatrix = new vector<Bucket<DATA> >[mutex_number];

	enumeration_started = false;
	cursed_iter = cursed_kmers.begin();
	current_minimizer = 0;
	current_overload=0;
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
	bucket.init_enumerator();
	while (bucket.has_next_kmer()) {
		bucket.next_kmer(kmer, small_minimizer);
		data=bucket.data_reserved_memory+kmer_id*sizeof(DATA);
		kint local_kmer=kmer.get_unhash_kmer_body(params.m,params.m_small,m_mask);
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

	// Cursed kmers
	if (kmer.multi_mini) {
		omp_set_lock(&multi_lock);
		cursed_kmers[kmer.kmer_s] = DATA();
		auto el=&(cursed_kmers[kmer.kmer_s]);
		omp_unset_lock(&multi_lock);
		return el;
	}

	// Transform the super minimizer to the used minimizer
	uint64_t small_minimizer = (((kmer.minimizer & mini_reduc_mask))>>(params.m-params.m_small));
	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Works because m - reduc <= 16
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	// omp_set_lock(&MutexBucket[mutex_idx]);
	uint32_t idx = bucket_indexes[matrix_idx];

	// GROGRO bucket check
	if (bucketMatrix[mutex_idx][idx-1].cleared) {
		omp_set_lock(&lock_overload[small_minimizer%bucket_overload]);
		kint local_kmer=kmer.get_unhash_kmer_body(params.m,params.m_small,m_mask);
		overload_kmers[small_minimizer%bucket_overload][local_kmer] = DATA();
		DATA* el=&overload_kmers[small_minimizer%bucket_overload][local_kmer];
		bucketMatrix[mutex_idx][idx-1].nb_kmers += 1;
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
	if ((not bucketMatrix[mutex_idx][idx-1].cleared) and bucketMatrix[mutex_idx][idx-1].nb_kmers>100 and bucketMatrix[mutex_idx][idx-1].nb_kmers*bucket_number/(total_number_kmers) >= GROGRO_THREASHOLD) {
		this->bucket_to_map(small_minimizer);
	}

	return value;
}



template <class DATA>
void DenseMenuYo<DATA>::insert_kmer_vector(const vector<kmer_full> & kmers,vector<bool>& newly_inserted) {
	vector<DATA *> result = this->get_kmer_vector(kmers);
	for(uint i(0);i<kmers.size(); ++i){
		if(result[i]==NULL){
			newly_inserted.push_back(true);
			bool newly_inserted = false;
			insert_kmer_no_mutex(kmers[i],newly_inserted,true);
		}else{
			newly_inserted.push_back(false);
		}
	}
	if (kmers[0].multi_mini) {
		return;
	}

	uint64_t small_minimizer = (((kmers[0].minimizer & mini_reduc_mask))>>(params.m-params.m_small));
	uint32_t mutex_idx = get_mutex(small_minimizer);
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);
	uint32_t idx = bucket_indexes[matrix_idx];

	// Bucket now too big to hold
	if ((not bucketMatrix[mutex_idx][idx-1].cleared) and bucketMatrix[mutex_idx][idx-1].nb_kmers>100 and bucketMatrix[mutex_idx][idx-1].nb_kmers*active_buckets/(total_number_kmers) >= GROGRO_THREASHOLD) {
		this->bucket_to_map(small_minimizer);
	}
	return;
}



template <class DATA>
DATA * DenseMenuYo<DATA>::insert_kmer_no_mutex(const kmer_full & kmer,bool& newly_inserted, bool already_checked) {
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
		auto el=& cursed_kmers[kmer.kmer_s];
		omp_unset_lock(&multi_lock);

		return el;
	}
	// Transform the super minimizer to the used minimizer

	uint64_t small_minimizer = (((kmer.minimizer & mini_reduc_mask))>>(params.m-params.m_small));
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
		active_buckets++;
	}

	// Grogro kmer
	if (bucketMatrix[mutex_idx][idx-1].cleared) {
		omp_set_lock(&lock_overload[small_minimizer%bucket_overload]);
		kint local_kmer=kmer.get_unhash_kmer_body(params.m,params.m_small,m_mask);
		overload_kmers[small_minimizer%bucket_overload][local_kmer] = DATA();
		DATA* el=&overload_kmers[small_minimizer%bucket_overload][local_kmer];
		omp_unset_lock(&lock_overload[small_minimizer%bucket_overload]);
		return (el);
	}

	// Insert the kmer in the right bucket
	DATA * value = bucketMatrix[mutex_idx][idx-1].insert_kmer(kmer);
	total_number_kmers++;

	return value;
}



template <class DATA>
DATA * DenseMenuYo<DATA>::get_kmer(kmer_full & kmer) {
	// Cursed kmers
	if (kmer.multi_mini) {
		omp_set_lock(&multi_lock);
		if (cursed_kmers.count(kmer.kmer_s) != 0) {
			auto el =&(cursed_kmers[kmer.kmer_s]);
			omp_unset_lock(&multi_lock);
			return el;
		} else {
			omp_unset_lock(&multi_lock);
			return (DATA *)NULL;
		}
	}
    uint64_t small_minimizer = (((kmer.minimizer & mini_reduc_mask))>>(params.m-params.m_small));

	kmer.minimizer_idx+=(params.m-params.m_small)/2;
	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Normal kmer
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	// omp_set_lock(&MutexBucket[mutex_idx]);
	uint32_t idx = bucket_indexes[matrix_idx];
	// No bucket
	if (idx == 0) {
		return (DATA *)NULL;
	}

	// GROGRO bucket check
	if (bucketMatrix[mutex_idx][idx-1].cleared) {
		kint local_kmer=kmer.get_unhash_kmer_body(params.m,params.m_small,m_mask);
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
	//WE ASSUME HERE THAT ALL KMERS are equally cursed
	if (kmers[0].multi_mini) {
		omp_set_lock(&multi_lock);
		for(uint i(0);i<kmers.size();++i){
			if (cursed_kmers.count(kmers[i].kmer_s) != 0) {
				result[i] =&cursed_kmers[kmers[i].kmer_s];
			}
		}
		omp_unset_lock(&multi_lock);
		return result;
	}
	//WE ASSUME HERE THAT ALL KMER HAVE THE SAME MINIMIZER
    uint64_t small_minimizer = (((kmers[0].minimizer & mini_reduc_mask))>>(params.m-params.m_small));

	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Normal kmer
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	uint32_t idx = bucket_indexes[matrix_idx];
	// No bucket
	if (idx == 0) {
		return result;
	}

	if (bucketMatrix[mutex_idx][idx-1].cleared) {
		omp_set_lock(&lock_overload[small_minimizer%bucket_overload]);
		for(uint i(0);i<kmers.size();++i){
			kint local_kmer=kmers[i].get_unhash_kmer_body(params.m,params.m_small,m_mask);
			if (overload_kmers[small_minimizer%bucket_overload].count(local_kmer) != 0) {
				result[i] =&overload_kmers[small_minimizer%bucket_overload][local_kmer];
			}else{
			}
			
		}
		omp_unset_lock(&lock_overload[small_minimizer%bucket_overload]);
		return result;
	}
	// Looks into the bucket for the right kmer
	result = bucketMatrix[mutex_idx][idx-1].find_kmer_vector(kmers);
	return result;
}



template <class DATA>
DATA * DenseMenuYo<DATA>::get_kmer_no_mutex(const kmer_full & kmer) {
	
	if (kmer.multi_mini) {
		omp_set_lock(&multi_lock);
		if (cursed_kmers.count(kmer.kmer_s) != 0) {
			auto el=&(cursed_kmers[kmer.kmer_s]);
			omp_unset_lock(&multi_lock);
			return el;
		} else {
			omp_unset_lock(&multi_lock);
			return (DATA *)NULL;
		}
	}
	
    uint64_t small_minimizer = (((kmer.minimizer & mini_reduc_mask))>>(params.m-params.m_small));

	uint32_t mutex_idx = get_mutex(small_minimizer);

	// Normal kmer
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	uint32_t idx = bucket_indexes[matrix_idx];
	// No bucket
	if (idx == 0) {
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
    uint64_t small_minimizer = (((kmer.minimizer & mini_reduc_mask))>>(params.m-params.m_small));
	uint32_t mutex_idx = get_mutex(small_minimizer);
	omp_set_lock(&MutexData[mutex_idx]);
}



template<class DATA>
void DenseMenuYo<DATA>::unprotect_data(const kmer_full & kmer) {
    uint64_t small_minimizer = (((kmer.minimizer & mini_reduc_mask))>>(params.m-params.m_small));
	uint32_t mutex_idx = get_mutex(small_minimizer);
	omp_unset_lock(&MutexData[mutex_idx]);
}



template<class DATA>
bool DenseMenuYo<DATA>::next(kmer_full & kmer) {
	if (not enumeration_started) {
		cursed_iter = cursed_kmers.begin();
		current_overload=0;
		overload_iter=overload_kmers[current_overload].begin();
		enumeration_started = true;
	}
	// Cursed kmers
	if (cursed_iter != cursed_kmers.end()) {
		kmer.kmer_s = cursed_iter->first;
		kmer.multi_mini = true;
		cursed_iter=std::next(cursed_iter);
		return true;
	}
	kmer.multi_mini =false;
	while(current_overload<bucket_overload){	
		if(overload_iter!=overload_kmers[current_overload].end()){
			kmer.kmer_s = overload_iter->first;
			
			bool reversed;
			kmer.minimizer=get_minimizer(kmer.kmer_s,params.k,kmer.minimizer_idx,params.m,reversed,kmer.multi_mini,params.mask_large_minimizer);
			kmer.hash_kmer_body(params.m,kmer.minimizer);
			overload_iter=std::next(overload_iter);
			return true;
		}
		current_overload++;
		if(current_overload<bucket_overload){
			overload_iter=overload_kmers[current_overload].begin();
		}
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

		// Already enumerated kmers
		if (idx != 0 and bucketMatrix[mutex_idx][idx-1].cleared) // Maybe a bug here
			idx = 0;
	}

	// Enumeration is over
	if (current_minimizer >= bucket_number)
		return false;

	if (not bucketMatrix[mutex_idx][idx-1].has_next_kmer()) {
		current_minimizer += 1;
		return this->next(kmer);
	}
	
	bucketMatrix[mutex_idx][idx-1].next_kmer(kmer, current_minimizer);
	kmer.minimizer_idx -= (params.m-params.m_small)/2;
	kmer.compute_mini(params.m);

	return true;
}



template<class DATA>
void DenseMenuYo<DATA>::restart_kmer_enumeration() {
	current_minimizer = 0;
}



ofstream bucket_size("bucketsize.txt");



static inline uint32_t mylog2(const uint32_t x) {
  uint32_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}



template<class DATA>
void DenseMenuYo<DATA>::stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & nb_cursed, uint64_t & largest_bucket) const {
	nb_buckets = 0;
	nb_kmers = 0;
	nb_skmers = 0;
	largest_bucket=0;
	bucket_size<<"size\n";

	//#pragma omp parallel for
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
			#pragma omp critical (filesize)
			{
				bucket_size<<(int)mylog2(bucketMatrix[mutex_idx][idx-1].skml.size())<<"\n";
			}
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
