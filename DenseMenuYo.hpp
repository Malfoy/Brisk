#include <iostream>
#include <cstdint>
#include <cstring>
#include <map>
#include <omp.h>
#include <vector>

#include <sys/resource.h>

#include "Kmers.hpp"
#include "buckets.hpp"

#ifndef DENSEMENUYO_H
#define DENSEMENUYO_H



template <class DATA>
class DenseMenuYo{
public:

	// Normal Buckets
	uint32_t * bucket_indexes;
	vector<Bucket<DATA> > * bucketMatrix;
	uint64_t k;
	// /!\ The number of bucket is the number of distinct minimizer / 4
	uint64_t minimizer_size;
	uint64_t bucket_number;
	uint64_t matrix_column_number;

	// Cursed kmer buckets
	robin_hood::unordered_map<kint, DATA> cursed_kmers;
	// Debug counts
	// robin_hood::unordered_flat_map<string, DATA> real_counts;


	DenseMenuYo(uint8_t k, uint8_t m, uint8_t minimizer_reduction=4);
	DATA * insert_kmer(kmer_full & kmer, const kint minimizer);
	DATA * get_kmer(kmer_full & kmer, const kint minimizer);
	// void add_kmers(vector<kmer_full>& v,uint64_t minimizer);
	// uint64_t dump_counting(bool check);
	// void dump_stats();

private:
	// Mutexes
	uint64_t mutex_order;
	uint64_t mutex_number;
	uint8_t m_reduc;
	// No idea on what that is ??
	// uint64_t call_ad;
	// uint64_t skm_total_size;
	// map<uint,uint> size_sk;

	// Number of mutex used
	vector<omp_lock_t> MutexWall;

	// Buffers
	kint buffered_minimizer;
	kmer_full * buffered_kmer;
	
	// uint64_t getMemorySelfMaxUsed();
};


template <class DATA>
DenseMenuYo<DATA>::DenseMenuYo(uint8_t k, uint8_t m, uint8_t minimizer_reduction){
	m_reduc = minimizer_reduction;
	// If m - reduc > 16, not enougth space for bucket idxs ! (because of uint32_t)
	assert(m - m_reduc <= 16);
	Bucket<DATA>::set_parameters(k, m-m_reduc);

	// Create the mutexes
	mutex_order = min((uint8_t)6,(uint8_t)(m-m_reduc));
	mutex_number = 1<<(2*mutex_order); /* DO NOT MODIFY, to increase/decrease modify mutex_order*/
	// Init the mutexes
	MutexWall.resize(mutex_number);
	for (uint64_t i(0); i < mutex_number; ++i) {
		omp_init_lock(&MutexWall[i]);
	}

	minimizer_size=m;
	this->k = k;
	bucket_number=1<<(2*(minimizer_size-m_reduc));
	matrix_column_number = bucket_number / mutex_number;
	bucket_indexes = new uint32_t[bucket_number];
	memset(bucket_indexes, 0, sizeof(uint32_t)*bucket_number);
	bucketMatrix = new vector<Bucket<DATA> >[mutex_number];

	buffered_minimizer = 0;
	buffered_kmer = NULL;
}


#define get_mutex(mini) (mini%mutex_number)
#define get_column(mini) (mini/mutex_number)
#define get_matrix_position(row_idx, col_idx) (row_idx * matrix_column_number + col_idx)


template <class DATA>
DATA * DenseMenuYo<DATA>::insert_kmer(kmer_full & kmer, const kint minimizer) {
	// if (buffered_kmer->kmer_s != kmer.kmer_s or minimizer != buffered_minimizer) {
		DATA * prev_val = this->get_kmer(kmer, minimizer);
		if (prev_val != NULL) {
			return prev_val;
		}
	// }

	// Cursed kmers
	if (kmer.multi_mini) {
		cursed_kmers[kmer.kmer_s] = DATA();
		return &(cursed_kmers[kmer.kmer_s]);
	}

	// Works because m - reduc <= 16
	uint32_t small_minimizer = (uint32_t)(minimizer >> (2 * m_reduc));
	uint32_t mutex_idx = get_mutex(small_minimizer);
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	omp_set_lock(&MutexWall[mutex_idx]);
	uint32_t idx = bucket_indexes[matrix_idx];
	// Create the bucket if not already existing
	if (idx == 0) {
		bucketMatrix[mutex_idx].emplace_back();
		bucket_indexes[matrix_idx] = bucketMatrix[mutex_idx].size();
		idx = bucket_indexes[matrix_idx];
	}
	// Insert the kmer in the right bucket
	DATA * value = bucketMatrix[mutex_idx][idx-1].insert_kmer(kmer);
	omp_unset_lock(&MutexWall[mutex_idx]);

	// buffered_minimizer = minimizer;
	// buffered_kmer = &kmer;
	return value;
}


template <class DATA>
DATA * DenseMenuYo<DATA>::get_kmer(kmer_full & kmer, const kint minimizer) {
	// Cursed kmers
	if (kmer.multi_mini) {
		if (cursed_kmers.count(kmer.kmer_s) != 0)
			return &(cursed_kmers[kmer.kmer_s]);
		else
			return (DATA *)NULL;
	}

	// Normal kmer
	uint32_t small_minimizer = (uint32_t)(minimizer >> (2 * m_reduc));
	uint32_t mutex_idx = get_mutex(small_minimizer);
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	omp_set_lock(&MutexWall[mutex_idx]);
	uint32_t idx = bucket_indexes[matrix_idx];
	// No bucket
	if (idx == 0) {
		omp_unset_lock(&MutexWall[mutex_idx]);
		return (DATA *)NULL;
	}

	// Looks into the bucket for the right kmer
	DATA * value = bucketMatrix[mutex_idx][idx-1].find_kmer(kmer);
	omp_unset_lock(&MutexWall[mutex_idx]);
	
	// buffered_minimizer = minimizer;
	// buffered_kmer = &kmer;

	return value;
}


#endif
