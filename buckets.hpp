#include <iostream>
#include <cstdint>
#include "Kmers.hpp"
#include "SuperKmerCount.hpp"



#ifndef BUCKETS_H
#define BUCKETS_H


//TODO CHANGE VECTOR of STruct TO Struct of vector
class Bucket{
public:
	vector<SKCL> skml;
	vector<uint8_t> values;
	uint32_t sorted_size;

	void insert_buffer();
	void add_kmers(vector<kmer_full>& kmers);
	void add_kmers_buffer( vector<kmer_full>& kmer);
	bool add_kmers_sorted( vector<kmer_full>& kmer);
	void print_kmers(string& result,const  string& mini)const ;
	uint64_t size()const;
	uint64_t number_kmer()const;
	Bucket(){sorted_size=0;}
};

#endif
