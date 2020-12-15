#include <iostream>
#include <cstdint>
#include <math.h>

#include "pow2.hpp"
#include "Kmers.hpp"
#include "parameters.hpp"



#ifndef SKL_H
#define SKL_H



class SKL {
public:
	uint32_t idx;
	uint32_t data_idx;
	uint8_t size;//1B
	uint8_t minimizer_idx;//1B
	/**
	 * The number of bytes that are really occupied by the nucleotides
	 */
	uint8_t bytes_used;


	SKL(kint kmer, const uint8_t mini_idx, uint32_t idx, uint8_t * nucleotides, uint32_t data_idx, const Parameters & params);
	SKL(const SKL & otto);
	SKL& operator=(const SKL& rhs);

	bool compact_right(const kmer_full & kmer, uint8_t * nucleotides, const Parameters & params);
	bool is_kmer_present(const kmer_full& kmer, uint8_t * nucleotides, const Parameters & params) const;

	// uint32_t interleaved_value(uint8_t * nucleotides, const Parameters & params)const;
	// uint32_t interleaved_value_max(uint8_t * nucleotides, uint max_fix_idx, const Parameters & params)const;

	uint suffix_size()const;
	uint prefix_size(const Parameters & params)const;
	kint get_prefix(const uint8_t * nucleotides, const Parameters & params)const;
	kint get_suffix(const uint8_t * nucleotides, const Parameters & params)const;

	void get_kmer(const uint8_t kmer_idx, const uint8_t * nucleotides, const kint & mini, kmer_full & kmer, const Parameters & params)const;
	// kint get_skmer(const uint8_t * nucleotides, const kint & mini)const;
	// kint get_compacted(const uint8_t * nucleotides)const;

	bool inf (const uint8_t * my_nucleotides, const SKL & skmer, const uint8_t * sk_nucleotides, const Parameters & params) const;
	int8_t interleaved_nucleotide(const uint nucl_idx, const uint8_t * nucleotides, const Parameters & params) const;

	void print(const uint8_t * nucleotides, const kint & mini, const Parameters & params) const;

private:
	kint get_right_overlap(uint8_t * nucleotides, const Parameters & params)const;
	// uint32_t param_interleaved_value(uint8_t * nucleotides, uint8_t baseval, const Parameters & params)const;

	kint get_compacted_kmer(const uint8_t kmer_idx, const uint8_t * nucleotides, const Parameters & params)const;
	/**
  * Get the nucleotide value at the position in parameter.
  * Position 0 is the first nucleotide of the prefix.
  * The minimizer nucleotides doesn't count.
  */
	uint byte_index(uint position, const Parameters & params)const;
};

/** Construct a superkmer from one kmer and the minimizer position.
	* @param kmer The unsigned int used to represent the binary kmer. The minimizer is not present
	* @param mini_idx The minimizer position in the kmer (equivalent to suffix size).
	*/
SKL::SKL(kint kmer, const uint8_t mini_idx, uint32_t idx, uint8_t * nucleotides, uint32_t data_idx, const Parameters & params) {
	this->idx = idx;
	this->data_idx = data_idx;
	// memset(nucleotides, 0, params.allocated_bytes);
	Pow2<kint> anc(2 * params.compacted_size - 8);

	for(uint i(0) ; i<(params.compacted_size/4) ; i++){
		nucleotides[params.allocated_bytes-1-i] = kmer / anc;
		kmer%=anc;
		anc>>=8;
	}

	if(params.compacted_size % 4 != 0){
		nucleotides[params.allocated_bytes-1-(params.compacted_size/4)] = (kmer << (2 * (4 - params.compacted_size % 4)));
	}

	// this->interleaved = 0;
	this->size = 1;
	this->minimizer_idx = mini_idx;
	// this->kmer_data = vector<DATA>(1);
	this->bytes_used=ceil(static_cast<float>(params.k - params.m_small)/4.);
};

SKL::SKL(const SKL & otto) {
	this->idx = otto.idx;
	// this->interleaved = otto.interleaved;
	this->size = otto.size;
	this->minimizer_idx = otto.minimizer_idx;
	this->data_idx = otto.data_idx;
	this->bytes_used = otto.bytes_used;
}

SKL& SKL::operator=(const SKL& rhs) {
	this->idx = rhs.idx;
	// this->interleaved = rhs.interleaved;
	this->size = rhs.size;
	this->minimizer_idx = rhs.minimizer_idx;
	this->data_idx = rhs.data_idx;
	this->bytes_used = rhs.bytes_used;

	return *this;
}

bool SKL::compact_right(const kmer_full & kmer, uint8_t * nucleotides, const Parameters & params) {
	if ((params.compacted_size + size) / 8 > params.allocated_bytes)
		return false;
	
	size_t kmer_suffix_size = kmer.minimizer_idx;
	size_t skm_suffix_size = this->suffix_size();
	// Suffix sizes are not matching
	if (kmer_suffix_size != skm_suffix_size+1) {
		return false;
	}
	
	// Get the rightest nucleotides of the skmer
	kint super_kmer_overlap = this->get_right_overlap(nucleotides, params);
	// Get the rightest nucleotides of the kmer to compact
	kint kmer_overlap = kmer.get_compacted(params.m_small);
	int nuc = kmer_overlap % 4;
	kmer_overlap >>= 2;
	kmer_overlap %= ((kint)1 << (2 * (params.k - 1 - params.m_small)));
	
	if(super_kmer_overlap == kmer_overlap){
		int byte_to_update = byte_index(params.compacted_size + size - 1, params);
		int padding = (4 - ((params.compacted_size + size) % 4)) % 4;
		nucleotides[byte_to_update] += (nuc << (2 * padding));
		size++;
		// Size of a kmer - size of minimizer + 1 for each supplementary nucleotide
		bytes_used = ceil(static_cast<float>(params.k - params.m_small + size - 1)/4.);
		minimizer_idx++;

		return true;
	}

	return false;
}


uint SKL::suffix_size() const{
	return this->minimizer_idx;
}


uint SKL::prefix_size(const Parameters & params) const{
	return (this->size + params.compacted_size - 1 - this->minimizer_idx);
}

kint SKL::get_prefix(const uint8_t * nucleotides, const Parameters & params) const {
	kint result = 0;
	if (prefix_size(params) == 0)
		return result;

	// Copy the compacted value
	int start = byte_index(prefix_size(params) - 1, params);
	int length_to_read = byte_index(0, params) - start + 1;
	memcpy(&result, &nucleotides[start], length_to_read);
	// Align the compacted value
	int offset = (4 - (prefix_size(params) % 4)) % 4;
	result >>= 2 * offset;

	return result;
}


kint SKL::get_suffix(const uint8_t * nucleotides, const Parameters & params)const {
	kint result = 0;
	if (suffix_size() == 0)
		return result;

	// Copy the compacted value
	int start = byte_index(prefix_size(params) + suffix_size() - 1, params);
	int length_to_read = byte_index(prefix_size(params), params) - start + 1;
	memcpy(&result, &nucleotides[start], length_to_read);
	// Align the compacted value
	int offset = (4 - ((prefix_size(params) + suffix_size()) % 4)) % 4;
	result >>= 2 * offset;

	return result;
}



kint SKL::get_right_overlap(uint8_t * nucleotides, const Parameters & params) const {
	kint result(this->get_compacted_kmer(size-1, nucleotides, params));

	result %= (kint)1 << (2 * (params.k - 1 - params.m_small));

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
uint SKL::byte_index(uint nucl_position, const Parameters & params) const{
	uint relative_pos = nucl_position / 4;
	return params.allocated_bytes - relative_pos - 1;
}


bool SKL::is_kmer_present(const kmer_full& kmer, uint8_t * nucleotides, const Parameters & params) const{
	if (kmer.minimizer_idx <= this->minimizer_idx and // Suffix long enougth
			kmer.minimizer_idx - this->minimizer_idx + size > 0) { // Prefix long enougth
		int kmer_idx = size - (this->minimizer_idx - kmer.minimizer_idx) - 1;
		return get_compacted_kmer(kmer_idx, nucleotides, params) == kmer.get_compacted(params.m_small);
	} else
		return false;
}



kint SKL::get_compacted_kmer(const uint8_t kmer_idx, const uint8_t * nucleotides, const Parameters & params) const {
	kint result = 0;

	// Copy the compacted value
	int start = byte_index(kmer_idx + params.compacted_size - 1, params);
	int length_to_read = byte_index(kmer_idx, params) - start + 1;
	memcpy(&result, &nucleotides[start], length_to_read);
	// Align the compacted value
	int offset = (4 - ((params.compacted_size + kmer_idx) % 4)) % 4;
	result = result >> (2 * offset);
	result = result & (((kint)1 << (2 * params.compacted_size)) - 1);

	return result;
}


void SKL::get_kmer(const uint8_t kmer_idx, const uint8_t * nucleotides, const kint & mini, kmer_full & kmer, const Parameters & params) const {
	kint compacted = get_compacted_kmer(kmer_idx, nucleotides, params);
	
	// Suffix preparation
	uint8_t suffix_size = this->minimizer_idx - (this->size - kmer_idx - 1);
	kint mask = ((kint)1 << (2 * suffix_size)) - 1;
	kint suffix = compacted & mask;

	// Prefix preparation
	mask = ~mask;
	kint prefix = (compacted & mask) << (2 * params.m_small);

	// Minimizer
	mask = ((kint)1 << (2 * params.m_small)) - 1;
	kmer.kmer_s = (mini & mask) << (2 * suffix_size);

	// Assemble everything
	kmer.kmer_s += prefix + suffix;
	kmer.multi_mini = false;
	kmer.minimizer_idx = suffix_size;
	kmer.minimizer = mini;
}


bool SKL::inf (const uint8_t * my_nucleotides, const SKL & skmer, const uint8_t * sk_nucleotides, const Parameters & params) const {
	// Compute my max number of iteration
	uint8_t my_suff_size = this->suffix_size();
	uint8_t my_pref_size = this->prefix_size(params);
	kint my_pref = this->get_prefix(my_nucleotides, params);
	kint my_suff = this->get_suffix(my_nucleotides, params);

	// Compute other max number of iteration
	uint8_t sk_suff_size = skmer.suffix_size();
	uint8_t sk_pref_size = skmer.prefix_size(params);
	kint sk_pref = skmer.get_prefix(sk_nucleotides, params);
	kint sk_suff = skmer.get_suffix(sk_nucleotides, params);

	// Compare when inside of both superkmers
	uint max_nucl = 2 * max(max(my_suff_size, my_pref_size), max(sk_suff_size, sk_pref_size));
	for (uint idx=0 ; idx<max_nucl ; idx++) {
		uint side_idx = idx / 2;
		uint8_t my_nucl = 0;
		uint8_t sk_nucl = 0;

		if (idx % 2 == 0) {
			// Verify borders
			if (side_idx >= my_suff_size and side_idx >= sk_suff_size)
				continue;
			else if (side_idx >= my_suff_size)
				return true;
			else if (side_idx >= sk_suff_size)
				return false;

			// Get the nucleotides to compare
			my_nucl = (my_suff >> (2 * (my_suff_size - 1 - side_idx))) & 0b11;
			sk_nucl = (sk_suff >> (2 * (sk_suff_size - 1 - side_idx))) & 0b11;
		} else {
			// Verify borders
			if (side_idx >= my_pref_size and side_idx >= sk_pref_size)
				continue;
			else if (side_idx >= my_pref_size)
				return true;
			else if (side_idx >= sk_pref_size)
				return false;

			// Get the nucleotides to compare
			my_nucl = (my_pref >> (2 * side_idx)) & 0b11;
			sk_nucl = (sk_pref >> (2 * side_idx)) & 0b11;
		}
		
		// Compare nucleotides
		if (my_nucl < sk_nucl)
			return true;
		else if (my_nucl > sk_nucl)
			return false;
	}
	return false;
}


int8_t SKL::interleaved_nucleotide(const uint nucl_idx, const uint8_t * nucleotides, const Parameters & params) const {
	uint side_idx = nucl_idx / 2;
	uint skmer_nucl_idx = 0; // From the beginning of the prefix

	if (nucl_idx % 2 == 0) {
		// Verify borders
		if (side_idx >= this->suffix_size())
			return -1;

		skmer_nucl_idx = this->prefix_size(params) + side_idx;
	} else {
		// Verify borders
		if (side_idx >= this->prefix_size(params))
			return -1;

		skmer_nucl_idx = this->prefix_size(params) - 1 - side_idx;
	}

	// Get the right bytes
	uint byte = nucleotides[byte_index(skmer_nucl_idx, params)];
	uint8_t val = (byte >> (2 * (3 - (skmer_nucl_idx % 4))))	 & 0b11;

	return val;
}



void SKL::print(const uint8_t * nucleotides, const kint & mini, const Parameters & params) const {
	kint prefix = get_prefix(nucleotides, params);
	kint suffix = get_suffix(nucleotides, params);

	print_kmer(prefix, prefix_size(params));
	print_kmer(mini, params.m_small);
	print_kmer(suffix, suffix_size());
}


#endif
