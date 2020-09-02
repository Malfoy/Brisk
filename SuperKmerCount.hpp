#include <iostream>
#include <cstdint>
#include "Kmers.hpp"



#ifndef SKC_H
#define SKC_H



class SKC {
public:
	skint sk;//8B or 16B
	uint32_t indice_value;//4B
	uint8_t size;//1B
	uint8_t minimizer_idx;//1B


	//TOTAL 16+4+2

	/** Construct a superkmer from one kmer and the minimizer position.
	* @param kmer The unsigned int 64 used to represent the binary kmer.
	* @param mini_idx The minimizer position in the kmer.
	*/
	SKC(const kint kmer, const uint8_t mini_idx,uint32_t indice_v);
	/** Look for the minimizer between the SKC and both of the fwd and rev kmer. Call the right compact function if needed.*/
	// bool add_kmer_old(const kmer_full& kmer);
	// bool add_kmer(const kmer_full& kmer);
	bool operator< (const  SKC& str);
	friend std::ostream & operator << (std::ostream& out, const SKC& skc);
	bool suffix_is_prefix(const kmer_full&);
	bool compact_right(const kint);
	bool compact_left(const kint);
	bool is_present(kint, uint64_t);
	bool is_present_brutforce(kmer_full kmer, uint8_t & mini_k_idx);
	bool compact_right(const kmer_full& kmf);
	bool is_present(kmer_full kmf);
	skint get_prefix()  const;
	skint get_suffix() const;
	bool query_kmer_bool(const kmer_full& kmer);
	uint32_t query_kmer_hash(const kmer_full& kmer);

	bool  operator < (const  SKC& str) const {
		kint s1((get_suffix()));
		kint s2((str.get_suffix()));
		if(minimizer_idx>=str.minimizer_idx){
			s1>>=(2*(minimizer_idx-str.minimizer_idx));
		}else{
			s2>>=(2*(str.minimizer_idx-minimizer_idx));
		}
		if(s1==s2){
			return minimizer_idx < str.minimizer_idx;
		}
		return  s1<s2;
	}
} __attribute__((packed));



#endif
