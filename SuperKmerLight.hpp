#include <iostream>
#include <cstdint>
#include <math.h>

#include "pow2.hpp"
#include "Kmers.hpp"



#ifndef SKCL_H
#define SKCL_H



template <class DATA>
class SKCL {
public:
	static void set_parameters(const uint8_t k, const uint8_t m);
	/**
	  * Array of bytes storing the superkmer.
	  * The storage is done using little endian encoding.
	  * The last byte store the first nucleotides of the superkmer.
	  * The array is filled backward when new kmers are compacted.
	  */
	static uint64_t allocated_bytes;
	uint32_t idx;
	uint32_t data_idx;
	uint32_t interleaved;
	uint8_t size;//1B
	uint8_t minimizer_idx;//1B
	/**
	 * The number of bytes that are really occupied by the nucleotides
	 */
	uint8_t bytes_used;


	SKCL(kint kmer, const uint8_t mini_idx, uint32_t idx, uint8_t * nucleotides, uint32_t data_idx);
	SKCL(const SKCL<DATA> & otto);
	SKCL& operator=(const SKCL& rhs);

	DATA * compact_right(const kmer_full & kmer, uint8_t * nucleotides, DATA * data);
	DATA * query_kmer(const kmer_full& kmer, uint8_t * nucleotides, DATA * data);

	uint32_t interleaved_value(uint8_t * nucleotides)const ;
	uint32_t interleaved_value_max(uint8_t * nucleotides, uint max_fix_idx=7)const;

	uint suffix_size()const;
	uint prefix_size()const;
	kint get_prefix(const uint8_t * nucleotides)const;
	kint get_suffix(const uint8_t * nucleotides)const;

	void get_kmer(const uint8_t kmer_idx, const uint8_t * nucleotides, const kint & mini, kmer_full & kmer)const;
	// kint get_skmer(const uint8_t * nucleotides, const kint & mini)const;
	// kint get_compacted(const uint8_t * nucleotides)const;

	void print(const uint8_t * nucleotides, const kint & mini) const;

	static uint8_t minimizer_size;
private:
	static uint8_t k;
	static uint8_t compacted_size;
	static kint compact_mask;
	
	kint get_right_overlap(uint8_t * nucleotides)const;
	uint32_t param_interleaved_value(uint8_t * nucleotides, uint8_t baseval)const;

	kint get_compacted_kmer(const uint8_t kmer_idx, const uint8_t * nucleotides)const;
	/**
  * Get the nucleotide value at the position in parameter.
  * Position 0 is the first nucleotide of the prefix.
  * The minimizer nucleotides doesn't count.
  */
	uint byte_index(uint position)const;
};


template <class DATA>
uint8_t SKCL<DATA>::k = 0;
template <class DATA>
uint8_t SKCL<DATA>::minimizer_size = 0;
template <class DATA>
uint64_t SKCL<DATA>::allocated_bytes= 0;
template <class DATA>
uint8_t SKCL<DATA>::compacted_size = 0;
template <class DATA>
kint SKCL<DATA>::compact_mask = 0;


template <class DATA>
void SKCL<DATA>::set_parameters(const uint8_t k, const uint8_t m) {
	SKCL::k = k;
	SKCL::minimizer_size = m;
	// SKCL::allocated_bytes=ceil(((float)(2*k-m))/4.);
	SKCL<DATA>::allocated_bytes = ceil((float)(2 * (k - m)) / 4.);
	SKCL::compacted_size = k - m;
	SKCL::compact_mask = (((kint)1) << (2*compacted_size)) - 1;
}



/** Construct a superkmer from one kmer and the minimizer position.
	* @param kmer The unsigned int used to represent the binary kmer. The minimizer is not present
	* @param mini_idx The minimizer position in the kmer (equivalent to suffix size).
	*/
template <class DATA>
SKCL<DATA>::SKCL(kint kmer, const uint8_t mini_idx, uint32_t idx, uint8_t * nucleotides, uint32_t data_idx) {
	this->idx = idx;
	this->data_idx = data_idx;
	memset(nucleotides, 0, SKCL<DATA>::allocated_bytes);
	Pow2<kint> anc(2 * SKCL<DATA>::compacted_size - 8);

	for(uint i(0) ; i<(SKCL<DATA>::compacted_size/4) ; i++){
		nucleotides[SKCL::allocated_bytes-1-i] = kmer / anc;
		kmer%=anc;
		anc>>=8;
	}

	if(compacted_size % 4 != 0){
		nucleotides[SKCL::allocated_bytes-1-(compacted_size/4)] = (kmer << (2 * (4 - compacted_size % 4)));
	}

	this->interleaved = 0;
	this->size = 1;
	this->minimizer_idx = mini_idx;
	// this->kmer_data = vector<DATA>(1);
	this->bytes_used=ceil(static_cast<float>(k - minimizer_size)/4.);
};

template <class DATA>
SKCL<DATA>::SKCL(const SKCL<DATA> & otto) {
	this->idx = otto.idx;
	this->interleaved = otto.interleaved;
	this->size = otto.size;
	this->minimizer_idx = otto.minimizer_idx;
	this->data_idx = otto.data_idx;
	this->bytes_used = otto.bytes_used;
}

template <class DATA>
SKCL<DATA>& SKCL<DATA>::operator=(const SKCL& rhs) {
	this->idx = rhs.idx;
	this->interleaved = rhs.interleaved;
	this->size = rhs.size;
	this->minimizer_idx = rhs.minimizer_idx;
	this->data_idx = rhs.data_idx;
	this->bytes_used = rhs.bytes_used;

	return *this;
}

template <class DATA>
DATA * SKCL<DATA>::compact_right(const kmer_full & kmer, uint8_t * nucleotides, DATA * data) {
	if ((compacted_size + size) / 8 > allocated_bytes)
		return NULL;
	
	size_t kmer_suffix_size = kmer.minimizer_idx;
	size_t skm_suffix_size = this->suffix_size();
	// Suffix sizes are not matching
	if (kmer_suffix_size != skm_suffix_size+1) {
		return NULL;
	}
	
	// Get the rightest nucleotides of the skmer
	kint super_kmer_overlap = this->get_right_overlap(nucleotides);
	// Get the rightest nucleotides of the kmer to compact
	kint kmer_overlap = kmer.get_compacted(SKCL<DATA>::minimizer_size);
	int nuc = kmer_overlap % 4;
	kmer_overlap >>= 2;
	kmer_overlap %= ((kint)1 << (2 * (k - 1 - minimizer_size)));
	
	if(super_kmer_overlap == kmer_overlap){
		int byte_to_update = SKCL<DATA>::byte_index(SKCL::compacted_size + size - 1);
		int padding = (4 - ((SKCL::compacted_size + size) % 4)) % 4;
		nucleotides[byte_to_update] += (nuc << (2 * padding));
		size++;
		// Size of a kmer - size of minimizer + 1 for each supplementary nucleotide
		bytes_used = ceil(static_cast<float>(k - minimizer_size + size - 1)/4.);
		minimizer_idx++;

		return data + this->size - 1;
	}

	return NULL;
}


template <class DATA>
uint SKCL<DATA>::suffix_size()const{
	return this->minimizer_idx;
}


template <class DATA>
uint SKCL<DATA>::prefix_size()const{
	return (this->size + SKCL::compacted_size - 1 - this->minimizer_idx);
}

template <class DATA>
kint SKCL<DATA>::get_prefix(const uint8_t * nucleotides)const {
	kint result = 0;
	if (prefix_size() == 0)
		return result;

	// Copy the compacted value
	int start = SKCL<DATA>::byte_index(prefix_size() - 1);
	int length_to_read = SKCL<DATA>::byte_index(0) - start + 1;
	memcpy(&result, &nucleotides[start], length_to_read);
	// Align the compacted value
	int offset = (4 - (prefix_size() % 4)) % 4;
	result >>= 2 * offset;

	return result;
}


template <class DATA>
kint SKCL<DATA>::get_suffix(const uint8_t * nucleotides)const {
	kint result = 0;
	if (suffix_size() == 0)
		return result;

	// Copy the compacted value
	int start = SKCL<DATA>::byte_index(prefix_size() + suffix_size() - 1);
	int length_to_read = SKCL<DATA>::byte_index(prefix_size()) - start + 1;
	memcpy(&result, &nucleotides[start], length_to_read);
	// Align the compacted value
	int offset = (4 - ((prefix_size() + suffix_size()) % 4)) % 4;
	result >>= 2 * offset;

	return result;
}



template <class DATA>
kint SKCL<DATA>::get_right_overlap(uint8_t * nucleotides)const {
	kint result(this->get_compacted_kmer(size-1, nucleotides));

	result %= (kint)1 << (2 * (k - 1 - minimizer_size));

	return result;
}


/**
  * Return the byte index corresponding to the nucletide position.
  * Position 0 is the first nucleotide of the prefix.
  * The minimizer nucleotides doesn't count.
  * 
  * @param nucl_position Nucleotide position in the sequence (0 for first prefix nucleotide)
  *
  * @return The Byte index in the datastructure.
  * The first 4 nucleotides are inside of the last byte of the byte array (little endian style).
  */
template <class DATA>
uint SKCL<DATA>::byte_index(uint nucl_position)const{
	// cout << "bytes " << (uint)allocated_bytes << " pos " << nucl_position << endl;
	// return (4 *SKCL<DATA>::allocated_bytes - nucl_position - 2) / 8;
	uint relative_pos = nucl_position / 4;
	return allocated_bytes - relative_pos - 1;
}


template <class DATA>
uint32_t SKCL<DATA>::param_interleaved_value(uint8_t * nucleotides, uint8_t baseval)const{
	baseval &= 0b11;
	kint prefix = get_prefix(nucleotides);
	kint suffix = get_suffix(nucleotides);
	uint8_t suf_size = suffix_size();
	uint8_t pref_size = prefix_size();

	uint32_t value = 0;
	// Suffix interleaved
	uint8_t max_suffix = min((uint)8, (uint)suf_size);
	uint8_t i;
	for (i=0 ; i<max_suffix ; i++) {
		// Get the value of the nucleotide at the position
		uint64_t nucl_value = (suffix >> (2 * (suf_size - i - 1))) & 0b11;
		// shift the value to the right place
		nucl_value <<= 30 - i*4;
		// Add the nucleotide to the interleaved
		value |= nucl_value;
	}
	// Fill missing suffix nucleotides
	for (; i<8 ; i++) {
		value |= ((uint64_t)baseval << (30 - i*4));
	}

 	// prefix interleaved
	uint8_t max_prefix = min((uint)8, (uint)pref_size);
	for (i=0 ; i<max_prefix ; i++) {
		// Get the value of the nucleotide at the position
		uint64_t nucl_value = (prefix >> (2 * i)) & 0b11;
		// shift the value to the right place
		nucl_value <<= 28 - i*4;
		// Add the nucleotide to the interleaved
		value |= nucl_value;
	}
	// Fill missing prefix nucleotides
	for (; i<8 ; i++) {
		value |= ((uint64_t)baseval << (28 - i*4));
	}

	return value;
}


template <class DATA>
uint32_t SKCL<DATA>::interleaved_value(uint8_t * nucleotides)const{
	return this->param_interleaved_value(nucleotides, 0);
}

template <class DATA>
uint32_t SKCL<DATA>::interleaved_value_max(uint8_t * nucleotides, uint max_fix_idx)const {
	uint32_t min_val = this->interleaved;
	
	uint prefix_size = min((uint)max_fix_idx, (uint)this->prefix_size());
	uint suffix_size = min((uint)max_fix_idx, (uint)this->suffix_size());
	uint fix_point;
	if (prefix_size <= suffix_size) {
		fix_point = min(max_fix_idx, 2 * (suffix_size - 1));
	} else {
		fix_point = min(max_fix_idx, 2 * (prefix_size - 1) + 1);
	}

	uint32_t increment_value = (uint32_t)1 << (30 - fix_point * 2);

	return min_val + increment_value;
}


template <class DATA>
DATA * SKCL<DATA>::query_kmer(const kmer_full& kmer, uint8_t * nucleotides, DATA * data) {
	int64_t start_idx  = (int64_t)this->minimizer_idx - (int64_t)kmer.minimizer_idx;
	if(start_idx<0 or (start_idx>=this->size)){
		return NULL;
	}

	uint64_t kmer_idx = size - start_idx - 1;
	if (get_compacted_kmer(kmer_idx, nucleotides) == kmer.get_compacted(minimizer_size)) {
		return data + kmer_idx;
	} else {
		return NULL;
	}
}


template <class DATA>
kint SKCL<DATA>::get_compacted_kmer(const uint8_t kmer_idx, const uint8_t * nucleotides)const {
	kint result = 0;

	// Copy the compacted value
	int start = SKCL<DATA>::byte_index(kmer_idx + compacted_size - 1);
	int length_to_read = SKCL<DATA>::byte_index(kmer_idx) - start + 1;
	memcpy(&result, &nucleotides[start], length_to_read);
	// Align the compacted value
	int offset = (4 - ((compacted_size + kmer_idx) % 4)) % 4;
	result >>= 2 * offset;
	result = result & (((kint)1 << (2 * compacted_size)) - 1);

	return result;
}


template <class DATA>
void SKCL<DATA>::get_kmer(const uint8_t kmer_idx, const uint8_t * nucleotides, const kint & mini, kmer_full & kmer)const {
	kint compacted = get_compacted_kmer(kmer_idx, nucleotides);
	
	// Suffix preparation
	uint8_t suffix_size = this->minimizer_idx - (this->size - kmer_idx - 1);
	kint mask = ((kint)1 << (2 * suffix_size)) - 1;
	kint suffix = compacted & mask;

	// Prefix preparation
	mask = ~mask;
	kint prefix = (compacted & mask) << (2 * SKCL<DATA>::minimizer_size);

	// Minimizer
	mask = ((kint)1 << (2 * SKCL<DATA>::minimizer_size)) - 1;
	kmer.kmer_s = (mini & mask) << (2 * suffix_size);

	// Assemble everything
	kmer.kmer_s += prefix + suffix;
	kmer.multi_mini = false;
	kmer.minimizer_idx = suffix_size;
}


template <class DATA>
void SKCL<DATA>::print(const uint8_t * nucleotides, const kint & mini) const {
	kint prefix = get_prefix(nucleotides);
	kint suffix = get_suffix(nucleotides);

	print_kmer(prefix, prefix_size());
	print_kmer(mini, minimizer_size);
	print_kmer(suffix, suffix_size());
}


#endif
