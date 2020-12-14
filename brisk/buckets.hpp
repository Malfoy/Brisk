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
	DATA * find_kmer_from_interleave(kmer_full& kmer, SKL & mockskm, uint8_t * mock_nucleotides);
	DATA * find_kmer_linear(const kmer_full& kmer, const uint64_t begin, const uint64_t end);
	DATA * find_kmer_log(const kmer_full & kmer, const uint64_t begin, const uint64_t end, const uint8_t nucleotide_idx);
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


template <class DATA>
void Bucket<DATA>::insert_buffer(){
	for(auto it(skml.begin()+sorted_size);it<skml.end();++it){
		it->interleaved=it->interleaved_value(
				this->nucleotides_reserved_memory + it->idx * params->allocated_bytes,
				*params
		);
	}

	sort(skml.begin()+sorted_size , skml.end(),[ ]( const SKL& lhs, const SKL& rhs ){
		return lhs.interleaved < rhs.interleaved;});
	inplace_merge(skml.begin(), skml.begin()+sorted_size ,skml.end() ,[ ]( const SKL& lhs, const SKL& rhs ){
		return lhs.interleaved < rhs.interleaved;});

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


static uint64_t max_i = 0;

template <class DATA>
DATA * Bucket<DATA>::find_kmer_from_interleave(kmer_full& kmer, SKL & mockskm, uint8_t * mock_nucleotides){
	DATA * data_pointer = NULL;

	uint64_t low = lower_bound(
		skml.begin(),
		skml.begin() + sorted_size,
		mockskm,
		[ ]( const SKL& lhs, const SKL& rhs ){
			return lhs.interleaved < rhs.interleaved;
		}
	) - skml.begin();
	uint32_t min_value = skml[low].interleaved;
	uint32_t value_max = mockskm.interleaved_value_max(mock_nucleotides, 5, *params);
	if (value_max < mockskm.interleaved)
		value_max = std::numeric_limits<uint32_t>::max();

	uint i = 0;
	while (data_pointer == NULL and low < (uint64_t)sorted_size and skml[low].interleaved < value_max) {
		i += 1;
		bool is_present = skml[low].is_kmer_present(kmer,
					this->nucleotides_reserved_memory + params->allocated_bytes * skml[low].idx,
					*params
		);
		if (is_present) {
			uint8_t kmer_position = skml[low].size - (skml[low].minimizer_idx - kmer.minimizer_idx) - 1;
			buffered_data = data_pointer;

			return this->data_reserved_memory + skml[low].data_idx + kmer_position;
		}

		low++;
	}
	if (i > max_i) {
		max_i = i;
		cout << "max linear " << i << " size " << skml.size() << endl;
		cout << min_value << endl;
		cout << value_max << endl;
		cout << (uint)kmer.minimizer_idx << " "; print_kmer(kmer.kmer_s, this->params->k); cout << endl;
	}

	return data_pointer;
}


template<class DATA>
DATA * Bucket<DATA>::find_kmer_log(const kmer_full & kmer, const uint64_t begin, const uint64_t end, const uint8_t nucleotide_idx) {
	// print_kmer(kmer.kmer_s, 63); cout << " " << begin << " " << end << " " << (uint)nucleotide_idx << endl;
	// Base case, linear interogation
	// If nucl idx 16 means end of the interleved value used to sort
	if (nucleotide_idx == 16 or end - begin < 1) {
		if (nucleotide_idx == 16) {
			uint8_t nucleotides_area[32]; // Max 2 * sizeof(kint)
			SKL mockskm(kmer.get_compacted(params->m_small), kmer.minimizer_idx, 0, nucleotides_area, 0, *params);
			mockskm.interleaved = mockskm.interleaved_value(nucleotides_area, *params);
			// #pragma omp critical
			// {
			// cout << "linear size " << (uint)(end - begin) << endl;
			// cout << mockskm.interleaved << " " << (uint)kmer.minimizer_idx << endl;
			// cout << skml[begin].interleaved << endl;
			// cout << skml[end].interleaved << endl;
			// }
		}
		return find_kmer_linear(kmer, begin, end);
	}

	// Select the middle skmer
	uint64_t middle = begin + (end - begin) / 2;
	SKL & skmer = skml[middle];
	uint8_t middle_nucl = skmer.interleaved >> (30 - (2 * nucleotide_idx));
	middle_nucl = middle_nucl & 0b11;

	// if exists, get the middle nucleotide of interest
	bool absent_nucl = true;
	auto ps_idx = nucleotide_idx / 2;
	uint8_t nucl = 0;
	if (nucleotide_idx % 2 == 0) { //Suffix nucleotide
		if (ps_idx < kmer.minimizer_idx) {
			nucl = kmer.kmer_s >> (2 * (kmer.minimizer_idx - 1 - ps_idx));
			nucl = nucl & 0b11;
			absent_nucl = false;
		}
	} else { // Prefix nucleotide
		if (ps_idx < this->params->k - this->params->m - kmer.minimizer_idx) {
			nucl = kmer.kmer_s >> (2 * (kmer.minimizer_idx + this->params->m_small + ps_idx));
			nucl = nucl & 0b11;
			absent_nucl = false;
		}
	}

	// Case 1 - Nucleotide is present
	if (not absent_nucl) {
		// cout << "case 1" << endl;
		if (nucl < middle_nucl) {
			if (begin < middle)
				return find_kmer_log(kmer, begin, middle-1, nucleotide_idx);
			else
				return NULL;
		} else if (nucl > middle_nucl) {
			if (middle < end)
				return find_kmer_log(kmer, middle+1, end, nucleotide_idx);
			else
				return NULL;
		} else { // equality
			uint8_t begin_nucl = (skml[begin].interleaved >> (30 - (2 * nucleotide_idx))) & 0b11;
			uint8_t end_nucl = (skml[end].interleaved >> (30 - (2 * nucleotide_idx))) & 0b11;

			// The current nucleotide is not determinant
			if (begin_nucl == end_nucl) {
				return find_kmer_log(kmer,begin, end, nucleotide_idx+1);
			}
			// Try the current skmer and recur both direction
			else {
				// debug_count += 1;
				DATA * ptr = find_kmer_linear(kmer, middle, middle);
				// This is the middle skmer !
				if (ptr != NULL)
					return ptr;
				// Recur left
				if (begin != middle)
					ptr = find_kmer_log(kmer, begin, middle-1, nucleotide_idx);
				if (ptr != NULL)
					return ptr;
				// Recur right
				if (middle != end)
					return find_kmer_log(kmer, middle+1, end, nucleotide_idx);
				return NULL;
			}
		}
	}

	// Case 2 - Can be any nucleotide
	else {
		uint8_t begin_nucl = (skml[begin].interleaved >> (30 - (2 * nucleotide_idx))) & 0b11;
		uint8_t end_nucl = (skml[end].interleaved >> (30 - (2 * nucleotide_idx))) & 0b11;
			// cout << "case 2" << endl;

		// The current nucleotide is not determinant
		if (begin_nucl == end_nucl) {
			return find_kmer_log(kmer,begin, end, nucleotide_idx+1);
		}
		// Try the current skmer and recur both direction
		else {
			// debug_count += 1;
			DATA * ptr = find_kmer_linear(kmer, middle, middle);
			// This is the middle skmer !
			if (ptr != NULL)
				return ptr;
			// Recur left
			if (middle != begin)
				ptr = find_kmer_log(kmer, begin, middle-1, nucleotide_idx);
			if (ptr != NULL)
				return ptr;
			// Recur right
			if (middle != end)
				return find_kmer_log(kmer, middle+1, end, nucleotide_idx);
			return NULL;
		}
	}

	cout << "WTF ??" << endl;
	return NULL;
}

template<class DATA>
DATA * Bucket<DATA>::find_kmer_linear(const kmer_full& kmer, const uint64_t begin, const uint64_t end) {
	// if (not (begin == end and begin == 0))
	// 	cout << begin << " " << end << endl;
	for (uint i=begin ; i<=end ; i++) {
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
	// for (uint i=sorted_size ; i<skml.size() ; i++) {
	// 	bool is_present = skml[i].is_kmer_present(
	// 			kmer,
	// 			this->nucleotides_reserved_memory + params->allocated_bytes * skml[i].idx,
	// 			*params
	// 	);

	// 	if (is_present) {
	// 		uint8_t kmer_position = skml[i].size - (skml[i].minimizer_idx - kmer.minimizer_idx) - 1;
	// 		buffered_data = this->data_reserved_memory + skml[i].data_idx + kmer_position;
	// 		return this->data_reserved_memory + skml[i].data_idx + kmer_position;
	// 	}
	// }

	// return NULL;
}


template <class DATA>
DATA * Bucket<DATA>::find_kmer(kmer_full& kmer) {
	debug_count = 0;
	if (sorted_size > 0) {
		DATA * ptr = find_kmer_log(kmer, 0, sorted_size-1, 0);

		// if (debug_count > 0) {
		// 	#pragma omp critical
		// 	{
		// 		cout << "dbg " << debug_count << endl;
		// 	}
		// }

		if (ptr != NULL)
			return ptr;
	}

	return find_kmer_unsorted(kmer);
}

// template <class DATA>
// DATA * Bucket<DATA>::find_kmer(kmer_full& kmer) {
// 	if (buffered_get == kmer.kmer_s) {
// 		return buffered_data;
// 	} else {
// 		buffered_get = kmer.kmer_s;
// 		buffered_data = NULL;
// 	}

// 	uint8_t nucleotides_area[32]; // Max 2 * sizeof(kint)
// 	SKL mockskm(kmer.get_compacted(params->m_small), kmer.minimizer_idx, 0, nucleotides_area, 0, *params);
// 	mockskm.interleaved = mockskm.interleaved_value(nucleotides_area, *params);

// 	// print_kmer(mockskm.interleaved >> 52, 6); cout << " interleaved" << endl;
	
// 	uint prefix_size = mockskm.prefix_size(*params);
// 	uint suffix_size = mockskm.suffix_size();

// 	uint size_interleave = min(prefix_size,suffix_size) * 2;
// 	// cout << "interleaved size " << size_interleave << endl;

// 	if(size_interleave>=6){
// 		DATA * val_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
// 		// cout << "INTERLEAVED OUT " << (uint64_t *)val_pointer << endl;
// 		if (val_pointer != NULL)
// 			return val_pointer;
// 	}else{
// 		if(suffix_size>prefix_size){
// 			//SUFFIX IS LARGER PREFIX IS MISSING
// 			if(size_interleave==4){
// 				for(uint64_t i(0);i<4;++i){
// 					mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 12);
// 					DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
// 					if(value_pointer != NULL){return value_pointer;}
// 					mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 12);
// 				}
// 			}
// 			if(size_interleave==2){
// 				for(uint64_t ii(0);ii<4;++ii){
// 					mockskm.interleaved+=ii<<(sizeof(mockskm.interleaved) * 8 - 8);
// 					for(uint64_t i(0);i<4;++i){
// 						mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 12);
// 						DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
// 						if(value_pointer != NULL){return value_pointer;}
// 						mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 12);
// 					}
// 					mockskm.interleaved-=ii<<(sizeof(mockskm.interleaved) * 8 - 8);
// 				}
// 			}
// 			if(size_interleave==0){
// 				//~ cout<<"0 prefix"<<endl;
// 				for(uint64_t iii(0);iii<4;++iii){
// 					mockskm.interleaved+=iii<<(sizeof(mockskm.interleaved) * 8 - 4);
// 					for(uint64_t ii(0);ii<4;++ii){
// 						mockskm.interleaved+=ii<<(sizeof(mockskm.interleaved) * 8 - 8);
// 						for(uint64_t i(0);i<4;++i){
// 							mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 12);
// 							DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
// 							if(value_pointer != NULL){return value_pointer;}
// 							mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 12);
// 						}
// 						mockskm.interleaved-=ii<<(sizeof(mockskm.interleaved) * 8 - 8);
// 					}
// 					mockskm.interleaved-=iii<<(sizeof(mockskm.interleaved) * 8 - 4);
// 				}
// 			}
// 		}else{
// 			//PREFIX IS LARGER, SUFFIX IS MISSING
// 			if(size_interleave==4){
// 				for(uint64_t i(0);i<4;++i){
// 					mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 10);
// 					DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
// 					if(value_pointer != NULL){return value_pointer;}
// 					mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 10);
// 				}
// 			}
// 			if(size_interleave==2){
// 				for(uint64_t ii(0);ii<4;++ii){
// 					mockskm.interleaved+=ii<<(sizeof(mockskm.interleaved) * 8 - 6);
// 					for(uint64_t i(0);i<4;++i){
// 						mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 10);
// 						DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
// 						if(value_pointer != NULL){return value_pointer;}
// 						mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 10);
// 					}
// 					mockskm.interleaved-=ii<<(sizeof(mockskm.interleaved) * 8 - 6);
// 				}
// 			}
// 			if(size_interleave==0){
// 				for(uint64_t iii(0);iii<4;++iii){
// 					mockskm.interleaved+=iii<<(sizeof(mockskm.interleaved) * 8 - 2);
// 					for(uint64_t ii(0);ii<4;++ii){
// 						mockskm.interleaved+=ii<<(sizeof(mockskm.interleaved) * 8 - 6);
// 						for(uint64_t i(0);i<4;++i){
// 							mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 10);
// 							DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
// 							if(value_pointer != NULL){return value_pointer;}
// 							mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 10);
// 						}
// 						mockskm.interleaved-=ii<<(sizeof(mockskm.interleaved) * 8 - 6);
// 					}
// 					mockskm.interleaved-=iii<<(sizeof(mockskm.interleaved) * 8 - 2);
// 				}
// 			}
// 		}
// 	}

// 	return find_kmer_unsorted(kmer);
// }


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
