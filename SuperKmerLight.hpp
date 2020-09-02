#include <iostream>
#include <cstdint>
#include "Kmers.hpp"



#ifndef SKCL_H
#define SKCL_H



class SKCL {
public:
	uint8_t nucleotides[byte_nuc+1];
	uint32_t indice_value;//4B
	uint8_t size;//1B
	uint8_t minimizer_idx;//1B
	uint8_t bytes_used;
	uint64_t interleaved;



	SKCL(kint kmer, const uint8_t mini_idx,uint32_t indice_v);
	bool  suffix_is_prefix(const SKCL& kmf)const;
	string get_string(const string& mini) const ;
	kint get_ith_kmer(uint ind)const;
	kint get_suffix()const;
	kint get_right_overlap()const;
	bool operator < (const  SKCL& str) const;
	bool query_kmer_bool(const kmer_full& kmer)const ;
	uint32_t query_kmer_hash(const kmer_full& kmer)const;
	bool compact_right(const kmer_full& kmf);
	bool suffix_is_prefix(const kmer_full& kmf)const;
	uint64_t interleaved_value()const;
	void print_all()const;
} __attribute__((packed));



#endif
