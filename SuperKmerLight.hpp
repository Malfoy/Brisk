#include <iostream>
#include <cstdint>
#include "Kmers.hpp"



#ifndef SKCL_H
#define SKCL_H



class SKCL {
private:
	static const uint64_t allocated_bytes=ceil(((float)(2*k-minimizer_size))/4.);
	/**
  * Return the byte index corresponding to the nucletide position.
  * Position 0 is the first nucleotide of the prefix.
  * The minimizer nucleotides doesn't count.
  * 
  * @param position Nucleotide position in the sequence
  *
  * @return The Byte index in the datastructure.
  * The first 4 nucleotides are inside of the last byte of the byte array (little endian style).
  */
	uint8_t get_nucleotide(uint8_t position)const;
	/**
  * Get the nucleotide value at the position in parameter.
  * Position 0 is the first nucleotide of the prefix.
  * The minimizer nucleotides doesn't count.
  */
	uint byte_index(uint position)const;

	static uint which_byte(uint i);
	uint nb_nucl()const;

public:
	/**
	  * Array of bytes storing the superkmer.
	  * The storage is done using little endian encoding.
	  * The last byte store the first nucleotides of the superkmer.
	  * The array is filled backward when new kmers are compacted.
	  */
	uint8_t nucleotides[SKCL::allocated_bytes];

	uint32_t indice_value;//4B
	uint8_t size;//1B
	uint8_t minimizer_idx;//1B
	/**
	 * The number of bytes that are really occupied by the nucleotides
	 */
	uint8_t bytes_used;
	uint64_t interleaved;



	SKCL(kint kmer, const uint8_t mini_idx,uint32_t indice_v);
	uint64_t interleaved_value()const ;
	uint64_t interleaved_value_max()const;
	bool  suffix_is_prefix(const SKCL& kmf)const;
	string get_string(const string& mini) const ;
	kint get_ith_kmer(uint ind)const;
	kint get_suffix()const;
	kint get_right_overlap()const;
	bool operator < (const  SKCL& str) const;
	bool query_kmer_bool(const kmer_full& kmer)const ;
	int32_t query_kmer_hash(const kmer_full& kmer)const;
	bool compact_right(const kmer_full& kmf);
	bool suffix_is_prefix(const kmer_full& kmf)const;
	void print_all()const;
	bool  is_lex_inferior(const SKCL& kmf)const;	
	uint suffix_size()const;
	uint prefix_size()const;
	uint interleaved_size()const;
} __attribute__((packed));



#endif
