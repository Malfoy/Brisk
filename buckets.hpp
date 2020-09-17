#include <iostream>
#include <cstdint>
#include <vector>

#include "SuperKmerLight.hpp"
#include "Kmers.hpp"



#ifndef BUCKETS_H
#define BUCKETS_H


//TODO CHANGE VECTOR of STruct TO Struct of vector
class Bucket{
private:
	static uint8_t k;
	static uint8_t minimizer_size;
public:
	static void set_parameters(uint8_t k, uint8_t m) {
		Bucket::k = k;
		Bucket::minimizer_size = m;
		SKCL::set_parameters(k, m);
	}


	vector<SKCL> skml;
	vector<uint8_t> values;
	uint32_t sorted_size;

	void insert_buffer();
	void add_kmers(vector<kmer_full>& kmers);
	void add_kmers_buffer( vector<kmer_full>& kmer);
	bool add_kmers_sorted( vector<kmer_full>& kmer);
	uint64_t print_kmers(string& result,const  string& mini, robin_hood::unordered_flat_map<string, uint8_t> & real_count, robin_hood::unordered_map<kint, uint8_t> * cursed_kmers, bool check)const;
	uint64_t size()const;
	uint64_t number_kmer()const;
	bool  find_kmer_from_interleave(kmer_full& kmer, SKCL& mockskm);
	bool find_kmer(kmer_full& kmer);
	uint64_t number_kmer_counted()const;
	Bucket(){sorted_size=0;}
};

#endif
