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
	// uint8_t * nucleotides; // TODO: GO EXTERN !
	uint32_t idx;
	uint32_t data_idx;
	uint8_t size;//1B
	uint8_t minimizer_idx;//1B
	/**
	 * The number of bytes that are really occupied by the nucleotides
	 */
	uint8_t bytes_used;
	uint64_t interleaved;

	// vector<DATA> kmer_data;


	SKCL(kint kmer, const uint8_t mini_idx, uint32_t idx, uint8_t * nucleotides, uint32_t data_idx);
	SKCL(const SKCL<DATA> & otto);
	SKCL& operator=(const SKCL& rhs) {return * this;};
	// ~SKCL();

	uint64_t interleaved_value(uint8_t * nucleotides)const ;
	uint64_t interleaved_value_max()const;
	// bool  suffix_is_prefix(const SKCL& kmf)const;
	// string get_string(const string& mini) const ;
	// kint get_suffix()const;
	// bool operator < (const  SKCL& str) const;
	// int32_t query_kmer_hash(const kmer_full& kmer)const;
	DATA * compact_right(const kmer_full & kmer, uint8_t * nucleotides, DATA * data);
	// bool suffix_is_prefix(const kmer_full& kmf)const;
	// void print_all()const;
	// bool  is_lex_inferior(const SKCL& kmf)const;	
	// uint interleaved_size()const;
	DATA * query_kmer(const kmer_full& kmer, uint8_t * nucleotides, DATA * data);
	uint suffix_size()const;
	uint prefix_size()const;

private:
	static uint8_t k;
	static uint8_t minimizer_size;
	static uint8_t compacted_size;
	static kint compact_mask;
	
	kint get_ith_kmer(uint idx, uint8_t * nucleotides)const;
	kint get_right_overlap(uint8_t * nucleotides)const;

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
	uint8_t get_nucleotide(uint8_t position, uint8_t * nucleotides)const;
	/**
  * Get the nucleotide value at the position in parameter.
  * Position 0 is the first nucleotide of the prefix.
  * The minimizer nucleotides doesn't count.
  */
	uint byte_index(uint position)const;

	static uint which_byte(uint i);
	uint nb_nucl()const;

}; //__attribute__((packed));


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
	SKCL::allocated_bytes=ceil(((float)(2*k-m))/4.);
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
	// this->nucleotides = otto.nucleotides;
	// memcpy(this->nucleotides, otto.nucleotides, SKCL<DATA>::allocated_bytes);

	this->idx = otto.idx;
	this->interleaved = otto.interleaved;
	this->size = otto.size;
	this->minimizer_idx = otto.minimizer_idx;
	// this->kmer_data = otto.kmer_data;
	this->data_idx = otto.data_idx;
	this->bytes_used = otto.bytes_used;
}

template <class DATA>
DATA * SKCL<DATA>::compact_right(const kmer_full & kmer, uint8_t * nucleotides, DATA * data) {
	size_t kmer_suffix_size = kmer.minimizer_idx;
	size_t skm_suffix_size = this->suffix_size();
	// Suffix sizes are not matching
	if (kmer_suffix_size != skm_suffix_size+1) {
		return NULL;
	}
	
	// Get the rightest nucleotides of the skmer
	kint super_kmer_overlap = this->get_right_overlap(nucleotides);
	// Get the rightest nucleotides of the kmer to compact
	kint kmer_overlap = kmer.get_compacted();
	int nuc = kmer_overlap % 4;
	kmer_overlap >>= 2;
	kmer_overlap %= ((kint)1 << (2 * (k - 1 - minimizer_size)));
	
	if(super_kmer_overlap == kmer_overlap){
		int byte_to_update = SKCL::allocated_bytes - 1 - (SKCL::compacted_size+size - 1) / 4;
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


// uint SKCL::interleaved_size()const{
// 	return min(prefix_size(),suffix_size())*2;
// }

template <class DATA>
uint SKCL<DATA>::which_byte(uint i){
	return SKCL<DATA>::allocated_bytes - ceil((float)i / 4);
}


template <class DATA>
kint SKCL<DATA>::get_ith_kmer(uint kmer_idx, uint8_t * nucleotides)const{
	kint result(0);
	
	int start = SKCL<DATA>::which_byte(SKCL::compacted_size + kmer_idx - 1);
	int length_to_read = SKCL<DATA>::which_byte(kmer_idx) - start + 1;
	
	memcpy(&result, &nucleotides[start], length_to_read);
	
	int offset = (4 - (compacted_size + kmer_idx - 1) % 4) % 4;
	result >>= 2 * offset;
	result &= compact_mask;
	return result;
}

template <class DATA>
kint SKCL<DATA>::get_right_overlap(uint8_t * nucleotides)const {
	//~ cout<<"right overlap"<<endl;
	kint result(this->get_ith_kmer(size, nucleotides));
	//~ print_kmer(result,31);cout<<endl;
	result %= (kint)1 << (2 * (k - 1 - minimizer_size));
	return result;
}



template <class DATA>
uint SKCL<DATA>::nb_nucl()const{
	return SKCL::compacted_size + this->size - 1;
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
	return SKCL<DATA>::allocated_bytes - 1 - nucl_position / 4;
}


/**
  * Get the nucleotide value at the position in parameter.
  * Position 0 is the first nucleotide of the prefix.
  * The minimizer nucleotides doesn't count.
  */
template <class DATA>
uint8_t SKCL<DATA>::get_nucleotide(uint8_t nucl_position, uint8_t * nucleotides)const {
	uint byte_pos = this->byte_index(nucl_position);
	//~ cout << " alloc " << SKCL::allocated_bytes << " B " << byte_pos<<" nucl position:	"<<(int)nucl_position<<endl;;
	//~ cout << "Byte: " << byte_pos << endl;
	uint8_t nucl = nucleotides[byte_pos];
	nucl >>= 2 * (3 - (nucl_position%4));
	nucl &= 0b11;
	return nucl;
}


template <class DATA>
uint64_t SKCL<DATA>::interleaved_value(uint8_t * nucleotides)const{
	uint64_t value = 0;
	// Suffix interleaved
	uint8_t max_suffix = min((uint)8, (uint)minimizer_idx);
	for (uint8_t i=0 ; i<max_suffix ; i++) {
		uint8_t nucl_position = this->nb_nucl() - minimizer_idx + i;
		// Get the value of the nucleotide at the position
		uint64_t nucl_value = this->get_nucleotide(nucl_position, nucleotides);
		// shift the value to the right place
		nucl_value <<= 62 - i*4;
		// Add the nucleotide to the interleaved
		value |= nucl_value;
	}

 	// prefix interleaved
	uint8_t max_prefix = min((uint)8,prefix_size());
	for (uint8_t i=0 ; i<max_prefix ; i++) {
		// Get the nucleotide position
		uint8_t nucl_position = this->nb_nucl() - minimizer_idx - i - 1;
		// Get the value of the nucleotide at the position
		uint64_t nucl_value = this->get_nucleotide(nucl_position, nucleotides);
		// shift the value to the right place
		nucl_value <<= 60 - i*4;
		// Add the nucleotide to the interleaved
		value |= nucl_value;
	}
	return value;
}

template <class DATA>
uint64_t SKCL<DATA>::interleaved_value_max()const {
	uint64_t value = interleaved;
	// Suffix interleaved
	for (uint64_t i=max(suffix_size(),(uint)3) ; i<8 ; i++) {
		value += ((uint64_t)3<<(62-4*i));
	}
	// prefix interleaved
	for (uint64_t i=max(prefix_size(),(uint)3) ; i<8 ; i++) {
		value += ((uint64_t)3<<(60-4*i));
	}
	return value;
}


template <class DATA>
DATA * SKCL<DATA>::query_kmer(const kmer_full& kmer, uint8_t * nucleotides, DATA * data) {
	// cout << "SKL - query_kmer" << endl;
	int64_t start_idx  = (int64_t)this->minimizer_idx - (int64_t)kmer.minimizer_idx;
	
	if(start_idx<0 or (start_idx>=this->size)){
		return NULL;
	}

	if (get_ith_kmer(size-start_idx, nucleotides) == kmer.get_compacted())//THE GET COMPACTED SHOULD BE MADE ABOVE
		return data + size - start_idx - 1;
	else
		return NULL;
}


// template <class DATA>
// int32_t SKCL<DATA>::query_kmer_hash(const kmer_full& kmer)const {
// 	//~ cout<<"query_kmer_hash"<<endl;	
// 	if(this->query_kmer_bool(kmer)){
// 		//~ cout<<"query kmer bool true"<<endl;
// 		return indice_value+kmer.minimizer_idx - (minimizer_idx-size+1);
// 	}
// 	//~ cout<<"query kmer bool FAIL"<<endl;
// 	return -1;
// }


#endif
