#include <iostream>
#include <cstdint>
#include <cstring>
#include <map>
#include <omp.h>
#include <vector>

#include <sys/resource.h>

// #include "Kmers.hpp"
#include "buckets.hpp"

#ifndef DENSEMENUYO_H
#define DENSEMENUYO_H



class DenseMenuYo{
public:
	// Mutexes
	uint64_t mutex_order;
	uint64_t mutex_number;

	// Normal Buckets
	uint32_t * bucket_indexes;
	vector<Bucket> * bucketMatrix;
	uint64_t k;
	uint64_t minimizer_size;
	uint64_t bucket_number;
	uint64_t matrix_column_number;

	// Cursed kmer buckets
	robin_hood::unordered_map<kint, uint8_t> cursed_kmers[1024];
	// Debug counts
	robin_hood::unordered_flat_map<string, uint8_t> real_counts;


	DenseMenuYo(uint8_t k, uint8_t m);
	void add_kmers(vector<kmer_full>& v,uint64_t minimizer);
	uint64_t dump_counting(bool check);
	void dump_stats();

private:
	// No idea on what that is ??
	uint64_t call_ad;
	uint64_t skm_total_size;
	// map<uint,uint> size_sk;

	// Number of mutex used
	vector<omp_lock_t> MutexWall;
	
	uint64_t getMemorySelfMaxUsed();
};


#endif
