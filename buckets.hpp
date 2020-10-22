#include <iostream>
#include <cstdint>
#include <vector>
#include <limits>

#include "SuperKmerLight.hpp"
#include "Kmers.hpp"



#ifndef BUCKETS_H
#define BUCKETS_H


//TODO CHANGE VECTOR of STruct TO Struct of vector
template <class DATA>
class Bucket{
public:
	static void set_parameters(uint8_t k, uint8_t m) {
		Bucket::k = k;
		Bucket::minimizer_size = m;
		SKCL::set_parameters(k, m);
	}

	vector<SKCL> skml;
	uint32_t sorted_size;

	uint nb_kmers;

	Bucket();
	Bucket(const Bucket<DATA> &) = delete;
	Bucket(Bucket<DATA> && bucket);
	~Bucket();

	Bucket<DATA> & operator=(const Bucket<DATA>&) = delete;
	Bucket<DATA> & operator=(Bucket<DATA>&& bucket);

	DATA * insert_kmer(kmer_full & kmer);
	DATA * find_kmer(kmer_full& kmer);

	void next_kmer(kmer_full & kmer, kint minimizer);
	bool has_next_kmer();

private:
	static uint8_t k;
	static uint8_t minimizer_size;

	SKCL * buffered_skmer;
	kint buffered_get;
	DATA * buffered_data;

	uint8_t * nucleotides_reserved_memory;
	uint32_t skmer_reserved;
	uint64_t next_data;
	uint64_t data_reserved_number;
	DATA * data_reserved_memory;

	uint32_t enumeration_skmer_idx;
	uint8_t enumeration_kmer_idx;

	DATA * find_kmer_unsorted(kmer_full& kmer);
	DATA * find_kmer_from_interleave(kmer_full& kmer, SKCL & mockskm, uint8_t * mock_nucleotides);
	DATA * insert_kmer_buffer(kmer_full & kmer);
	void insert_buffer();
	void data_space_update();
};

// Static variable definition
template <class DATA>
uint8_t Bucket<DATA>::k = 0;
template <class DATA>
uint8_t Bucket<DATA>::minimizer_size = 0;


template <class DATA>
Bucket<DATA>::Bucket() {
	this->sorted_size = 0;
	this->buffered_skmer = (SKCL *)NULL;
	this->buffered_get = (kint)1 << (sizeof(kint) * 8 - 1);
	this->buffered_data = 0;

	this->nb_kmers = 0;

	this->skml.reserve(10);
	this->nucleotides_reserved_memory = (uint8_t *)malloc(SKCL::allocated_bytes * skml.capacity());
	this->next_data = 0;
	this->data_reserved_number = 10;
	this->data_reserved_memory = (DATA *)malloc(sizeof(DATA) * this->data_reserved_number);

	this->enumeration_skmer_idx = 0;
	this->enumeration_kmer_idx = 0;
}

template <class DATA>
Bucket<DATA>::Bucket(Bucket<DATA> && bucket)
: skml( std::move(bucket.skml) )
, sorted_size( bucket.sorted_size )
, nb_kmers( bucket.nb_kmers )
, buffered_skmer( bucket.buffered_skmer )
, buffered_get( bucket.buffered_get )
, buffered_data( bucket.buffered_data )
, nucleotides_reserved_memory( bucket.nucleotides_reserved_memory )
, skmer_reserved( bucket.skmer_reserved )
, next_data( bucket.next_data )
, data_reserved_number( bucket.data_reserved_number )
, data_reserved_memory( bucket.data_reserved_memory )
, enumeration_skmer_idx( bucket.enumeration_skmer_idx )
, enumeration_kmer_idx( bucket.enumeration_kmer_idx )
{
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
void Bucket<DATA>::data_space_update() {
	if (this->next_data == this->data_reserved_number) {
		size_t to_reserve = (int)min(1000., 0.5 * this->data_reserved_number);
		this->data_reserved_memory = (DATA *)realloc(this->data_reserved_memory, sizeof(DATA) * (this->data_reserved_number + to_reserve));
		this->data_reserved_number += to_reserve;
		
		this->buffered_data = NULL;
		this->buffered_get = (kint)1 << (sizeof(kint) * 8 - 1);
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
						SKCL::allocated_bytes * buffered_skmer->idx
		);

		if (is_compacted) {
			this->nb_kmers += 1;
			this->next_data += 1;
			return data_reserved_memory + this->next_data - 1;
		}
	}
	
	// 2 - Sort if needed
	if(skml.size()-sorted_size>00){
		this->insert_buffer();
	}

	// 3 - Create a new skmer
	DATA * value = this->insert_kmer_buffer(kmer);
	this->nb_kmers += 1;
	this->next_data += 1;

	return value;
}


template <class DATA>
void Bucket<DATA>::insert_buffer(){
	for(auto it(skml.begin()+sorted_size);it<skml.end();++it){
		it->interleaved=it->interleaved_value(
				this->nucleotides_reserved_memory + it->idx * SKCL::allocated_bytes);
	}

	sort(skml.begin()+sorted_size , skml.end(),[ ]( const SKCL& lhs, const SKCL& rhs ){
		return lhs.interleaved < rhs.interleaved;});
	inplace_merge(skml.begin(), skml.begin()+sorted_size ,skml.end() ,[ ]( const SKCL& lhs, const SKCL& rhs ){
		return lhs.interleaved < rhs.interleaved;});

	buffered_skmer = NULL;
	sorted_size=skml.size();
}


template <class DATA>
DATA * Bucket<DATA>::insert_kmer_buffer(kmer_full & kmer){
	DATA * value_pointer = NULL;
	this->data_space_update();

	// Scale Superkmer vector capacity if needed
	if(skml.size()==skml.capacity()){
		skml.reserve(skml.capacity()*1.5);
		this->nucleotides_reserved_memory = (uint8_t *)realloc(
										this->nucleotides_reserved_memory,
										skml.capacity() * SKCL::allocated_bytes
		);
	}

	// Create a new superkmer
	skml.emplace_back(
		kmer.get_compacted(SKCL::minimizer_size),
		(int)kmer.minimizer_idx,
		skml.size(),
		this->nucleotides_reserved_memory + (this->skml.size() * SKCL::allocated_bytes),
		this->next_data
	);
	buffered_skmer = &(skml[skml.size()-1]);
	value_pointer = this->data_reserved_memory + this->next_data;

	return value_pointer;
}


template <class DATA>
DATA * Bucket<DATA>::find_kmer_from_interleave(kmer_full& kmer, SKCL & mockskm, uint8_t * mock_nucleotides){
	DATA * data_pointer = NULL;

	uint64_t low = lower_bound(
		skml.begin(),
		skml.begin() + sorted_size,
		mockskm,
		[ ]( const SKCL& lhs, const SKCL& rhs ){
			return lhs.interleaved < rhs.interleaved;
		}
	) - skml.begin();
	uint32_t value_max = mockskm.interleaved_value_max(mock_nucleotides, 5);
	if (value_max < mockskm.interleaved)
		value_max = std::numeric_limits<uint32_t>::max();
	
	while (data_pointer == NULL and low < (uint64_t)sorted_size and skml[low].interleaved < value_max) {
		int8_t kmer_position = skml[low].query_kmer(kmer,
					this->nucleotides_reserved_memory + SKCL::allocated_bytes * skml[low].idx
		);
		if (kmer_position >= 0) {
			data_pointer = this->data_reserved_memory + skml[low].data_idx + kmer_position;
		}

		low++;
	}

	buffered_data = data_pointer;
	return data_pointer;
}


template <class DATA>
DATA * Bucket<DATA>::find_kmer_unsorted(kmer_full& kmer) {
	int8_t kmer_position = -1;
	uint i = 0;

	for (i=sorted_size ; i<skml.size() ; i++) {
		kmer_position = skml[i].query_kmer(
				kmer,
				this->nucleotides_reserved_memory + SKCL::allocated_bytes * skml[i].idx
		);

		if (kmer_position >= 0)
			break;
	}

	if (kmer_position >= 0) {
		buffered_data = this->data_reserved_memory + skml[i].data_idx + kmer_position;
		return this->data_reserved_memory + skml[i].data_idx + kmer_position;
	} else
		return NULL;
}


template <class DATA>
DATA * Bucket<DATA>::find_kmer(kmer_full& kmer) {
	if (buffered_get == kmer.kmer_s) {
		return buffered_data;
	} else {
		buffered_get = kmer.kmer_s;
	}

	// print_kmer(kmer.kmer_s, 17); cout << endl;
	static uint8_t * nucleotides_area = new uint8_t[SKCL::allocated_bytes];
	SKCL mockskm(kmer.get_compacted(SKCL::minimizer_size), kmer.minimizer_idx, 0, nucleotides_area, 0);
	mockskm.interleaved = mockskm.interleaved_value(nucleotides_area);

	// print_kmer(mockskm.interleaved >> 52, 6); cout << " interleaved" << endl;
	
	uint prefix_size = mockskm.prefix_size();
	uint suffix_size = mockskm.suffix_size();

	uint size_interleave = min(prefix_size,suffix_size) * 2;
	// cout << "interleaved size " << size_interleave << endl;

	if(size_interleave>=6){
		DATA * val_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
		// cout << "INTERLEAVED OUT " << (uint64_t *)val_pointer << endl;
		if (val_pointer != NULL)
			return val_pointer;
	}else{
		if(suffix_size>prefix_size){
			//SUFFIX IS LARGER PREFIX IS MISSING
			if(size_interleave==4){
				for(uint64_t i(0);i<4;++i){
					mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 12);
					DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
					if(value_pointer != NULL){return value_pointer;}
					mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 12);
				}
			}
			if(size_interleave==2){
				for(uint64_t ii(0);ii<4;++ii){
					mockskm.interleaved+=ii<<(sizeof(mockskm.interleaved) * 8 - 8);
					for(uint64_t i(0);i<4;++i){
						mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 12);
						DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
						if(value_pointer != NULL){return value_pointer;}
						mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 12);
					}
					mockskm.interleaved-=ii<<(sizeof(mockskm.interleaved) * 8 - 8);
				}
			}
			if(size_interleave==0){
				//~ cout<<"0 prefix"<<endl;
				for(uint64_t iii(0);iii<4;++iii){
					mockskm.interleaved+=iii<<(sizeof(mockskm.interleaved) * 8 - 4);
					for(uint64_t ii(0);ii<4;++ii){
						mockskm.interleaved+=ii<<(sizeof(mockskm.interleaved) * 8 - 8);
						for(uint64_t i(0);i<4;++i){
							mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 12);
							DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
							if(value_pointer != NULL){return value_pointer;}
							mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 12);
						}
						mockskm.interleaved-=ii<<(sizeof(mockskm.interleaved) * 8 - 8);
					}
					mockskm.interleaved-=iii<<(sizeof(mockskm.interleaved) * 8 - 4);
				}
			}
		}else{
			//PREFIX IS LARGER, SUFFIX IS MISSING
			if(size_interleave==4){
				for(uint64_t i(0);i<4;++i){
					mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 10);
					DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
					if(value_pointer != NULL){return value_pointer;}
					mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 10);
				}
			}
			if(size_interleave==2){
				for(uint64_t ii(0);ii<4;++ii){
					mockskm.interleaved+=ii<<(sizeof(mockskm.interleaved) * 8 - 6);
					for(uint64_t i(0);i<4;++i){
						mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 10);
						DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
						if(value_pointer != NULL){return value_pointer;}
						mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 10);
					}
					mockskm.interleaved-=ii<<(sizeof(mockskm.interleaved) * 8 - 6);
				}
			}
			if(size_interleave==0){
				for(uint64_t iii(0);iii<4;++iii){
					mockskm.interleaved+=iii<<(sizeof(mockskm.interleaved) * 8 - 2);
					for(uint64_t ii(0);ii<4;++ii){
						mockskm.interleaved+=ii<<(sizeof(mockskm.interleaved) * 8 - 6);
						for(uint64_t i(0);i<4;++i){
							mockskm.interleaved+=i<<(sizeof(mockskm.interleaved) * 8 - 10);
							DATA * value_pointer = find_kmer_from_interleave(kmer, mockskm, nucleotides_area);
							if(value_pointer != NULL){return value_pointer;}
							mockskm.interleaved-=i<<(sizeof(mockskm.interleaved) * 8 - 10);
						}
						mockskm.interleaved-=ii<<(sizeof(mockskm.interleaved) * 8 - 6);
					}
					mockskm.interleaved-=iii<<(sizeof(mockskm.interleaved) * 8 - 2);
				}
			}
		}
	}

	return find_kmer_unsorted(kmer);
}


template <class DATA>
bool Bucket<DATA>::has_next_kmer() {
	if (enumeration_skmer_idx >= skml.size())
		return false;

	SKCL & skmer = skml[enumeration_skmer_idx];
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

	SKCL & skmer = skml[enumeration_skmer_idx];
	skmer.get_kmer(
			enumeration_kmer_idx,
			this->nucleotides_reserved_memory + (skmer.idx * SKCL::allocated_bytes),
			minimizer,
			kmer
	);
	
	enumeration_kmer_idx += 1;
}


#endif
