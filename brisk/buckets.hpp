#include <iostream>
#include <cstdint>
#include <vector>
#include <limits>

#include "SuperKmerLight.hpp"
#include "Kmers.hpp"



#ifndef BUCKETS_H
#define BUCKETS_H


template <class DATA>
class Bucket{
public:
	vector<SKL> skml;
	uint32_t sorted_size;

	uint nb_kmers;

	Bucket(Parameters * params);
	Bucket(const Bucket<DATA> &) = delete;
	Bucket(Bucket<DATA> && bucket);
	~Bucket();

	Bucket<DATA> & operator=(const Bucket<DATA>&) = delete;
	Bucket<DATA> & operator=(Bucket<DATA>&& bucket);

	DATA * insert_kmer(kmer_full & kmer);
	DATA * find_kmer(kmer_full& kmer);

	void next_kmer(kmer_full & kmer, kint minimizer);
	bool has_next_kmer();

	void print();

private:
	friend class BriskWriter;

	Parameters * params;

	SKL * buffered_skmer;
	kint buffered_get;
	DATA * buffered_data;

	uint8_t * nucleotides_reserved_memory;
	uint32_t skmer_reserved;
	uint64_t next_data;
	uint64_t data_reserved_number;
	DATA * data_reserved_memory;

	uint32_t enumeration_skmer_idx;
	uint8_t enumeration_kmer_idx;

	uint debug_count;

	bool debug;

	DATA * find_kmer_unsorted(kmer_full& kmer);
	// DATA * find_kmer_from_interleave(kmer_full& kmer, SKL & mockskm, uint8_t * mock_nucleotides);
	DATA * find_kmer_linear(const kmer_full& kmer, const int64_t begin, const int64_t end);
	DATA * find_kmer_log(const kmer_full & kmer, const int64_t begin, const int64_t end, const uint nucleotide_idx);
	DATA * insert_kmer_buffer(kmer_full & kmer);
	void insert_buffer();
	void data_space_update();
};


template <class DATA>
Bucket<DATA>::Bucket(Parameters * params) {
	this->sorted_size = 0;
	this->buffered_skmer = (SKL *)NULL;
	this->buffered_get = ((kint)1) << (sizeof(kint) * 8 - 1);
	this->buffered_data = 0;

	this->params = params;

	this->nb_kmers = 0;

	this->skml.reserve(10);
	this->nucleotides_reserved_memory = (uint8_t *)malloc(params->allocated_bytes * skml.capacity());
	memset(this->nucleotides_reserved_memory, 0, params->allocated_bytes * skml.capacity());
	this->next_data = 0;
	this->data_reserved_number = 10;
	this->data_reserved_memory = (DATA *)malloc(sizeof(DATA) * this->data_reserved_number);

	this->enumeration_skmer_idx = 0;
	this->enumeration_kmer_idx = 0;

	this->debug = false;
}

template <class DATA>
Bucket<DATA>::Bucket(Bucket<DATA> && bucket)
: skml( std::move(bucket.skml) )
, sorted_size( bucket.sorted_size )
, nb_kmers( bucket.nb_kmers )
, params( bucket.params )
, buffered_skmer( bucket.buffered_skmer )
, buffered_get( ((kint)1) << (sizeof(kint) * 8 - 1) )
, buffered_data( bucket.buffered_data )
, nucleotides_reserved_memory( bucket.nucleotides_reserved_memory )
, skmer_reserved( bucket.skmer_reserved )
, next_data( bucket.next_data )
, data_reserved_number( bucket.data_reserved_number )
, data_reserved_memory( bucket.data_reserved_memory )
, enumeration_skmer_idx( bucket.enumeration_skmer_idx )
, enumeration_kmer_idx( bucket.enumeration_kmer_idx )
, debug( bucket.debug )
{
	if (debug)
			cout << "move construct modified" << endl;
	bucket.nucleotides_reserved_memory = NULL;
	bucket.data_reserved_memory = NULL;
}

template <class DATA>
Bucket<DATA> & Bucket<DATA>::operator=(Bucket<DATA>&& bucket) {
	cout << "move assign !!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
}


template <class DATA>
Bucket<DATA>::~Bucket() {
	free(this->nucleotides_reserved_memory);
	free(this->data_reserved_memory);
}


template <class DATA>
void Bucket<DATA>::print() {
	for (uint i=0 ; i<sorted_size ; i++) {
		print_kmer(
			skml[i].get_prefix(nucleotides_reserved_memory + (skml[i].idx*params->allocated_bytes), *params),
			skml[i].prefix_size(*params)
		); cout << " ";
		print_kmer(
			skml[i].get_suffix(nucleotides_reserved_memory + (skml[i].idx*params->allocated_bytes), *params),
			skml[i].suffix_size()
		); cout << endl;
	}
}


template <class DATA>
void Bucket<DATA>::data_space_update() {
	if (this->next_data == this->data_reserved_number) {
		auto factor = this->data_reserved_number > 500 ? 0.2 : 0.5;
		size_t to_reserve = (int)(factor * this->data_reserved_number);
		this->data_reserved_memory = (DATA *)realloc(this->data_reserved_memory, sizeof(DATA) * (this->data_reserved_number + to_reserve));
		this->data_reserved_number += to_reserve;
		
		if (debug)
			cout << "realloc modified" << endl;
		this->buffered_data = NULL;
		this->buffered_get = ((kint)1) << (sizeof(kint) * 8 - 1);
	}
}


template <class DATA>
DATA * Bucket<DATA>::insert_kmer(kmer_full & kmer) {
	// 0 - Update space for DATA if needed
	this->data_space_update();
	// 1 - Try to compact with the last kmer
	if (buffered_skmer != NULL) {
		bool is_compacted = buffered_skmer->compact_right(
				kmer,
				this->nucleotides_reserved_memory +
						this->params->allocated_bytes * buffered_skmer->idx,
				*params
		);

		if (is_compacted) {
			this->nb_kmers += 1;
			this->next_data += 1;
			
			this->buffered_data = data_reserved_memory + this->next_data - 1;
			this->buffered_get = kmer.kmer_s;
			return buffered_data;
		}
	}
	
	// 2 - Sort if needed
	if(skml.size()-sorted_size>10){
		this->insert_buffer();
	}

	// 3 - Create a new skmer
	DATA * value = this->insert_kmer_buffer(kmer);
	this->nb_kmers += 1;
	this->next_data += 1;
	this->buffered_data = value;
	this->buffered_get = kmer.kmer_s;

	return value;
}

// bool inf (const uint8_t * my_nucleotides, const SKL & skmer, const uint8_t * sk_nucleotides, const Parameters & params)

template <class DATA>
void Bucket<DATA>::insert_buffer(){
	Parameters * params = this->params;
	uint8_t * nucleotides_reserved_memory = this->nucleotides_reserved_memory;

	auto comp_function = [&](const SKL & a, const SKL & b) {
		bool val = a.inf(
				nucleotides_reserved_memory + a.idx * params->allocated_bytes,
				b,
				nucleotides_reserved_memory + b.idx * params->allocated_bytes,
				*params
		);
		// cout << val << endl;
		return val;
	};

	sort(skml.begin()+sorted_size, skml.end(), comp_function);
	inplace_merge(skml.begin(), skml.begin()+sorted_size, skml.end(), comp_function);

	buffered_skmer = NULL;
	this->buffered_data = NULL;
	this->buffered_get = ((kint)1) << (sizeof(kint) * 8 - 1);
	sorted_size=skml.size();
}


template <class DATA>
DATA * Bucket<DATA>::insert_kmer_buffer(kmer_full & kmer){
	DATA * value_pointer = NULL;
	this->data_space_update();

	// Scale Superkmer vector capacity if needed
	if(skml.size()==skml.capacity()){
		auto old_capacity = skml.capacity();
		auto factor = old_capacity > 500 ? 1.2 : 1.5;
		skml.reserve(skml.capacity()*factor);

		this->nucleotides_reserved_memory = (uint8_t *)realloc(
				this->nucleotides_reserved_memory,
				skml.capacity() * params->allocated_bytes
		);

		auto diff_capacity = skml.capacity() - old_capacity;
		memset(
				this->nucleotides_reserved_memory + params->allocated_bytes * old_capacity,
				0,
				params->allocated_bytes * diff_capacity
		);
	}

	// Create a new superkmer
	skml.emplace_back(
		kmer.get_compacted(params->m_small),
		(int)kmer.minimizer_idx,
		skml.size(),
		this->nucleotides_reserved_memory + (this->skml.size() * params->allocated_bytes),
		this->next_data,
		*params
	);
	buffered_skmer = &(skml[skml.size()-1]);
	value_pointer = this->data_reserved_memory + this->next_data;

	return value_pointer;
}


template<class DATA>
DATA * Bucket<DATA>::find_kmer_log(const kmer_full & kmer, const int64_t begin, const int64_t end, const uint nucleotide_idx) {
	if (debug)
		cout << begin << " " << end << " " << nucleotide_idx << endl;
	// Base case, linear interogation
	if (end - begin < 1) {
		return find_kmer_linear(kmer, begin, end);
	}

	uint side_idx = nucleotide_idx / 2;
	uint max_side_idx = max(kmer.suffix_size(), kmer.prefix_size(params->k, params->m_small));
	if (side_idx > max_side_idx)
		return find_kmer_linear(kmer, begin, end);

	// Select the middle skmer
	uint64_t middle = begin + (end - begin) / 2;
	SKL & skmer = skml[middle];

	uint8_t nucl = 0;
	// Suffix work
	if (nucleotide_idx % 2 == 0) {
		// If outiside of kmer, no information on this nucleotide
		if (side_idx >= kmer.suffix_size()) {
			if (debug)
				cout << "suffix too short" << endl;
			bool same_beginning = true;
			for (uint i=0 ; i<= nucleotide_idx ; i++) {
				int8_t begin_nucl = skml[begin].interleaved_nucleotide(i,
							this->nucleotides_reserved_memory + params->allocated_bytes * skml[begin].idx,
							*params);
				int8_t end_nucl = skml[end].interleaved_nucleotide(i,
							this->nucleotides_reserved_memory + params->allocated_bytes * skml[end].idx,
							*params);

				if (begin_nucl != end_nucl) {
					same_beginning = false;
					break;
				}
			}

			if (same_beginning)
				return this->find_kmer_log(kmer, begin, end, nucleotide_idx + 1);
			else {
				DATA * val = this->find_kmer_log(kmer, begin, middle, 0);
				if (val == NULL)
					val = this->find_kmer_log(kmer, middle + 1, end, 0);
				return val;
			}
		}
		// Must be after the middle
		else if (side_idx >= skmer.suffix_size())
			return this->find_kmer_log(kmer, middle+1, end, 0);
		// Get the nucleotide of interest of the kmer
		nucl = (kmer.kmer_s >> (2 * (kmer.minimizer_idx - 1 - side_idx))) & 0b11;

	// Prefix work
	} else {
		// If outiside of kmer, no information on this nucleotide
		if (side_idx >= kmer.prefix_size(params->k, params->m_small)) {
			if (debug)
				cout << "prefix too short" << endl;
			bool same_beginning = true;
			for (uint i=0 ; i<= nucleotide_idx ; i++) {
				int8_t begin_nucl = skml[begin].interleaved_nucleotide(i,
							this->nucleotides_reserved_memory + params->allocated_bytes * skml[begin].idx,
							*params);
				int8_t end_nucl = skml[end].interleaved_nucleotide(i,
							this->nucleotides_reserved_memory + params->allocated_bytes * skml[end].idx,
							*params);

				if (begin_nucl != end_nucl) {
					same_beginning = false;
					break;
				}
			}

			if (same_beginning)
				return this->find_kmer_log(kmer, begin, end, nucleotide_idx + 1);
			else {
				DATA * val = this->find_kmer_log(kmer, begin, middle, 0);
				if (val == NULL)
					val = this->find_kmer_log(kmer, middle + 1, end, 0);
				return val;
			}
		}
		// Must be after the middle
		else if (side_idx >= skmer.prefix_size(*params))
			return this->find_kmer_log(kmer, middle+1, end, 0);
		// Get the nucleotide of interest of the kmer
		nucl = (kmer.kmer_s >> (2 * (kmer.minimizer_idx + this->params->m_small + side_idx))) & 0b11;
	}
	if (debug)
		cout << "nucl " << (uint)nucl << endl;

	// Get the nucleotide of interest for the middle superkmer of the list
	uint8_t middle_nucl = skmer.interleaved_nucleotide(
				nucleotide_idx,
				this->nucleotides_reserved_memory + params->allocated_bytes * skmer.idx,
				*params
	);
	if (debug)
		cout << "m nucl " << (uint)middle_nucl << endl;

	if (nucl < middle_nucl) {
		return find_kmer_log(kmer, begin, middle-1, 0);
	} else if (nucl > middle_nucl) {
		return find_kmer_log(kmer, middle+1, end, 0);
	} else { // equality
		return find_kmer_log(kmer,begin, end, nucleotide_idx+1);
	}

	cout << "WTF ??" << endl;
	return NULL;
}

template<class DATA>
DATA * Bucket<DATA>::find_kmer_linear(const kmer_full& kmer, const int64_t begin, const int64_t end) {
	for (int i=begin ; i<=end ; i++) {
		bool is_present = skml[i].is_kmer_present(
				kmer,
				this->nucleotides_reserved_memory + params->allocated_bytes * skml[i].idx,
				*params
		);
		debug_count += 1;

		if (is_present) {
			uint8_t kmer_position = skml[i].size - (skml[i].minimizer_idx - kmer.minimizer_idx) - 1;
			buffered_data = this->data_reserved_memory + skml[i].data_idx + kmer_position;
			return this->data_reserved_memory + skml[i].data_idx + kmer_position;
		}
	}

	return NULL;
}


template <class DATA>
DATA * Bucket<DATA>::find_kmer_unsorted(kmer_full& kmer) {
	return find_kmer_linear(kmer, sorted_size, skml.size()-1);
}


template <class DATA>
DATA * Bucket<DATA>::find_kmer(kmer_full& kmer) {
	if (sorted_size > 0) {
		DATA * ptr = find_kmer_log(kmer, 0, sorted_size-1, 0);

		if (ptr != NULL)
			return ptr;
	}

	return find_kmer_unsorted(kmer);
}


template <class DATA>
bool Bucket<DATA>::has_next_kmer() {
	if (enumeration_skmer_idx >= skml.size())
		return false;

	SKL & skmer = skml[enumeration_skmer_idx];
	if (enumeration_kmer_idx >= skmer.size) {	
		enumeration_skmer_idx += 1;
		enumeration_kmer_idx = 0;
		return has_next_kmer();
	}	

	return true;
}


template <class DATA>
void Bucket<DATA>::next_kmer(kmer_full & kmer, kint minimizer) {
	// Nothing to do here
	if (not has_next_kmer())
		return;

	SKL & skmer = skml[enumeration_skmer_idx];
	skmer.get_kmer(
			enumeration_kmer_idx,
			this->nucleotides_reserved_memory + (skmer.idx * params->allocated_bytes),
			minimizer,
			kmer,
			*params
	);
	
	enumeration_kmer_idx += 1;
}


#endif
