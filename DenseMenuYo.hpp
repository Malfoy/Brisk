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
	vector<Bucket> * bucketMatrix;
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
	
	// uint64_t getMemorySelfMaxUsed();
};


template <class DATA>
DenseMenuYo<DATA>::DenseMenuYo(uint8_t k, uint8_t m, uint8_t minimizer_reduction){
	m_reduc = minimizer_reduction;
	// If m - reduc > 16, not enougth space for bucket idxs ! (because of uint32_t)
	assert(m - m_reduc <= 16);
	Bucket::set_parameters(k, m-m_reduc);

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
	bucketMatrix = new vector<Bucket>[mutex_number];

}


#define get_mutex(mini) (mini%mutex_number)
#define get_column(mini) (mini/mutex_number)
#define get_matrix_position(row_idx, col_idx) (row_idx * matrix_column_number + col_idx)


template <class DATA>
DATA * DenseMenuYo<DATA>::insert_kmer(kmer_full & kmer, const kint minimizer) {
	// Works because m - reduc <= 16
	uint32_t small_minimizer = (uint32_t)(minimizer >> (2 * m_reduc));
	uint32_t mutex_idx = get_mutex(small_minimizer);
	uint32_t column_idx = get_column(small_minimizer);
	uint64_t matrix_idx = get_matrix_position(mutex_idx, column_idx);

	omp_set_lock(&MutexWall[mutex_idx]);
	uint32_t idx = bucket_indexes[matrix_idx];
	if (idx == 0) {
		bucketMatrix[mutex_idx].push_back(Bucket());
		bucket_indexes[matrix_idx] = bucketMatrix[mutex_idx].size();
		idx = bucket_indexes[matrix_idx];
	}
	// Quick working code TODO: replace
	vector<kmer_full> v;
	v.push_back(kmer);
	bucketMatrix[mutex_idx][idx-1].add_kmers(v);
	// /TODO
	omp_unset_lock(&MutexWall[mutex_idx]);

	cout << "TODO: DenseMenuYo - insert_kmer" << endl;
	return NULL;
}


#endif
