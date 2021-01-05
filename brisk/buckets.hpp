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
	DATA * find_kmer_linear(kmer_full& kmer, const int64_t begin, const int64_t end);
	DATA * find_kmer_log(kmer_full & kmer);
	DATA * find_kmer_log(kmer_full & kmer, const int64_t begin, const int64_t end, const uint8_t nucleotide_idx);
	inline DATA * find_recur_log_split(kmer_full & kmer, const int64_t begin, const int64_t end, const uint8_t nucleotide_idx);
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
	if(skml.size()-sorted_size>2){
		// cout << "Sorting" << endl;
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
inline DATA * Bucket<DATA>::find_recur_log_split(kmer_full & kmer, const int64_t begin, const int64_t end, const uint8_t nucleotide_idx) {
	int8_t begin_nucl = skml[begin].interleaved_nucleotide(nucleotide_idx,
				this->nucleotides_reserved_memory + params->allocated_bytes * skml[begin].idx,
				*params);
	int8_t end_nucl = skml[end].interleaved_nucleotide(nucleotide_idx,
				this->nucleotides_reserved_memory + params->allocated_bytes * skml[end].idx,
				*params);
	// if (debug)
	// 	cout << "start/end " << (int)begin_nucl << " " << (int)end_nucl << endl;

	uint middle = begin + (end - begin) / 2;
	if (begin_nucl == end_nucl) {
		return find_kmer_log(kmer, begin, end, nucleotide_idx + 1);
	} else {
		DATA * val = find_kmer_log(kmer, begin, middle, nucleotide_idx);
		if (val != NULL)
			return val;
		return find_kmer_log(kmer, middle+1, end, nucleotide_idx);
	}
}


template<class DATA>
DATA * Bucket<DATA>::find_kmer_log(kmer_full & kmer, const int64_t begin, const int64_t end, const uint8_t nucleotide_idx) {
	// if (debug)
	// 	cout << begin << " " << end << " " << (uint)nucleotide_idx << endl;
	// Base case, linear interogation
	if (end - begin < 1) {
		return find_kmer_linear(kmer, begin, end);
	}

	uint8_t side_idx = nucleotide_idx / 2;
	uint8_t max_side_idx = max(kmer.suffix_size(), kmer.prefix_size(params->k, params->m_small));
	if (side_idx > max_side_idx or nucleotide_idx == 255) {
		// cout << "find max " << begin << " " << end << endl;
		return find_kmer_linear(kmer, begin, end);
	}

	int8_t nucl = kmer.interleaved_nucleotide(nucleotide_idx, params->k, params->m_small, debug);
	// if (debug)
	// 	cout << "nucl " << (int)nucl << endl;
	
	if (nucl < 0) {
		return find_recur_log_split(kmer, begin, end, nucleotide_idx);
	} else {
		// Select the middle skmer
		uint64_t middle = begin + (end - begin) / 2;
		SKL & mid_skmer = skml[middle];
		int8_t middle_nucl = mid_skmer.interleaved_nucleotide(nucleotide_idx,
							this->nucleotides_reserved_memory + params->allocated_bytes * mid_skmer.idx,
							*params);
		// if (debug) {
		// 	cout << "mid " << (int)middle_nucl << endl << endl;
		// }

		// Middle nucleotide absent
		if (middle_nucl < 0) {
			return find_kmer_log(kmer, middle+1, end, nucleotide_idx);
		} else {
			if (nucl < middle_nucl) {
				return find_kmer_log(kmer, begin, middle-1, nucleotide_idx);
			} else if (nucl > middle_nucl) {
				return find_kmer_log(kmer, middle+1, end, nucleotide_idx);
			} else {
				return find_recur_log_split(kmer, begin, end, nucleotide_idx);
			}
		}
	}

}

class find_params {
public:
	int64_t begin;
	int64_t end;
	int8_t start_letter;
	int8_t stop_letter;
	uint8_t start_interleaved_idx;
	uint8_t current_interleaved_idx;

	find_params(uint64_t b, uint64_t e, int8_t sl, int8_t spl, uint8_t ii, uint8_t cii) :
		begin(b), end(e), start_letter(sl), stop_letter(spl), start_interleaved_idx(ii), current_interleaved_idx(cii)
	{};

	void print() {
		cout << begin << " " << end << " " << (int)start_letter << " " << (int)stop_letter << " " << (uint)start_interleaved_idx << " " << (uint)current_interleaved_idx << endl;
	}
};


template<class DATA>
DATA * Bucket<DATA>::find_kmer_log(kmer_full & kmer) {
	// DEBUG PRINT
	// for (SKL & skmer : skml) {
	// 	skmer.print(nucleotides_reserved_memory + params->allocated_bytes * skmer.idx, (kint)1, *params); cout << endl;
	// }
	// cout << endl;
	// Prepare the kmer to search
	// cout << "find log" << endl;
	// kmer.print(params->k, params->m_small);

	vector<int8_t> kmer_interleved = kmer.compute_interleaved(params->k, params->m_small);
	// for (auto val : kmer_interleved) {
	// 	cout << (int)val << "\t";
	// }
	// cout << endl;

	vector<int8_t> current_interleved(kmer_interleved);


	// Prepare the search parameters
	int begin = 0; int end = sorted_size - 1;
	vector<int> heap;
	uint interleaved_size = kmer_interleved.size();
	vector<int8_t> middle_interleved(interleaved_size, -2);

	// Log search
	while (begin <= end or heap.size() > 0) {
		// cout << "new loop " << begin << "\t" << end << endl;
		// Base case - not found, restore previous search from heap
		if (begin > end) {
			// cout << "Restore ctx" << endl;
			uint interleaved_idx = heap.back(); heap.pop_back();
			end = heap.back(); heap.pop_back();
			begin = heap.back(); heap.pop_back();

			// First time modified
			if (current_interleved[interleaved_idx] == -2) {
				// cout << "All -1" << endl;
				// Set all unknown values to -1
				for (uint idx=interleaved_idx ; idx<interleaved_size ; idx += 2)
					current_interleved[idx] = -1;
			}
			// Restore the unknown following -2 and increment the current interleave
			else if (current_interleved[interleaved_idx] == -1) {
				// cout << "restore -2" << endl;
				current_interleved[interleaved_idx] = 0;
				for (uint idx=interleaved_idx+2 ; idx<interleaved_size ; idx += 2)
					current_interleved[idx] = -2;
			}
			// In other cases, change the value of the interleaved
			else {
				// cout << "increase current" << endl;
				current_interleved[interleaved_idx] += 1;
				for (uint idx=interleaved_idx+2 ; idx<interleaved_size ; idx += 2)
					if (kmer_interleved[idx] == -2)
						current_interleved[idx] = -2;
			}

			// If remaining possible values, save the context
			if (current_interleved[interleaved_idx] < 3) {
				heap.push_back(begin); heap.push_back(end); heap.push_back(interleaved_idx);
			}

			// cout << "boundaries " << begin << "\t" << end << endl;
			// cout << "current interleaved" << endl;
			// for (auto val : current_interleved) {
			// 	cout << (int)val << "\t";
			// }
			// cout << endl;
			// cout << "Heap" << endl;
			// for (uint idx=0 ; idx<heap.size() ; idx += 3)
			// 	cout << heap[idx] << "\t" << heap[idx+1] << "\t" << heap[idx+2] << endl;
			// cout << endl;
		}

		// Log case - get the middle superkmer
		int middle = begin + (end - begin) / 2;
		SKL & mid_skmer = skml[middle];
		mid_skmer.compute_interleaved(middle_interleved, nucleotides_reserved_memory + params->allocated_bytes * mid_skmer.idx, *params);

		// cout << "middle:" << endl;
		// mid_skmer.print(nucleotides_reserved_memory + params->allocated_bytes * mid_skmer.idx, (kint)1, *params); cout << endl;
		// for (auto val : middle_interleved) {
		// 	cout << (int)val << "\t";
		// }
		// cout << endl << endl;

		// Compare middle superkmer and searched superkmer
		bool found = true;
		for (uint idx=0 ; found and idx<interleaved_size ; idx++) {
			// Undefined current nucleotide
			if (current_interleved[idx] == -2) {
				// cout << "Unknown interleved idx " << idx << endl;
				// Save context for restoration
				heap.push_back(begin); heap.push_back(end); heap.push_back((int)idx);
				// Nothing to do here
				end = begin - 1;
				found = false;
			}
			// Before middle
			else if (current_interleved[idx] < middle_interleved[idx]) {
				// cout << "before middle " << idx << endl;
				end = middle - 1;
				found = false;
			}
			// After middle
			else if (current_interleved[idx] > middle_interleved[idx]) {
				// cout << "after middle " << idx << endl;
				begin = middle + 1;
				found = false;
			}
		}
		// cout << endl;

		if (found) {
			// cout << "FOUND §§" << endl;
			uint kmer_position = mid_skmer.prefix_size(*params) - kmer.prefix_size(params->k, params->m_small);
			// cout << (uint *)(this->data_reserved_memory + mid_skmer.data_idx + kmer_position) << endl;
			return this->data_reserved_memory + mid_skmer.data_idx + kmer_position;
		}

		// cout << "Heap" << endl;
		// for (uint idx=0 ; idx<heap.size() ; idx += 3)
		// 	cout << heap[idx] << "\t" << heap[idx+1] << "\t" << heap[idx+2] << endl;
		// cout << endl;
	}

	return NULL;
}

// template<class DATA>
// DATA * Bucket<DATA>::find_kmer_log(kmer_full & kmer) {
// 	cout << endl << "Find log" << endl;
// 	vector<find_params> call_stack;

// 	for (SKL & skmer : skml) {
// 		skmer.print(nucleotides_reserved_memory + params->allocated_bytes * skmer.idx, (kint)1, *params); cout << endl;
// 	}
// 	cout << endl;

// 	// Init the call stack
// 	int8_t first_nucl = kmer.interleaved_nucleotide(0, params->k, params->m_small, debug);	
// 	call_stack.emplace_back(0, sorted_size-1, first_nucl, first_nucl == -1 ? 3 : first_nucl, 0, 0);

// 	while (not call_stack.empty()) {
// 		// Get the current parameters
// 		find_params boundaries = call_stack[call_stack.size()-1];
// 		call_stack.pop_back();
// 		boundaries.print();

// 		// Base case
// 		if (boundaries.end - boundaries.begin < 1 or boundaries.current_interleaved_idx == 255) {
// 			DATA * val = find_kmer_linear(kmer, boundaries.begin, boundaries.end);
// 			if (val != NULL)
// 				return val;
// 		}
// 		else {
// 			// Get the letter of the middle superkmer
// 			uint64_t middle_idx = boundaries.begin + (boundaries.end - boundaries.begin) / 2;
// 			SKL & mid_skmer = skml[middle_idx];
// 			int8_t mid_letter = mid_skmer.interleaved_nucleotide(boundaries.current_interleaved_idx,
// 								nucleotides_reserved_memory + params->allocated_bytes * mid_skmer.idx,
// 								*params);
// 			cout << "mid skmer " << middle_idx << " " << (int) mid_letter << endl;

// 			// Check the boundaries
// 			if (boundaries.start_letter <= mid_letter and mid_letter <= boundaries.stop_letter) {
// 				// Multiple possible letter and try to divide into multiple requests
// 				if (boundaries.start_letter != boundaries.stop_letter) {
// 					// First letters before mid
// 					for (int8_t letter=boundaries.start_letter ; letter<mid_letter ; letter++) {
// 						call_stack.emplace_back(
// 							boundaries.begin, middle_idx-1,
// 							letter, letter,
// 							boundaries.start_interleaved_idx, boundaries.start_interleaved_idx);
// 					}
// 					// Last letters after mid
// 					for (int8_t letter=mid_letter+1 ; letter<=boundaries.stop_letter ; letter++) {
// 						call_stack.emplace_back(
// 							middle_idx+1, boundaries.end,
// 							letter, letter,
// 							boundaries.start_interleaved_idx, boundaries.start_interleaved_idx);
// 					}
// 					// Letter of mid
// 					call_stack.emplace_back(
// 						boundaries.begin, boundaries.end,
// 						mid_letter, mid_letter,
// 						boundaries.start_interleaved_idx, boundaries.current_interleaved_idx);
// 				}
// 				// Single possible letter
// 				else {
// 					while (mid_letter == boundaries.start_letter and mid_letter == boundaries.stop_letter) {
// 						boundaries.start_letter = kmer.interleaved_nucleotide(
// 							boundaries.current_interleaved_idx+1, params->k, params->m, debug);
// 						boundaries.stop_letter =
// 					}
// 				}
// 			}
// 			// Before mid
// 			else if (boundaries.stop_letter < mid_letter) {

// 			}
// 			// After mid
// 			else /*if (mid_letter < boundaries.start_letter)*/ {

// 			}
// 		}

// 		cout << "Stack" << endl;
// 		for (auto & p : call_stack) {
// 			p.print();
// 		}
// 		exit(0);
// 	}

// 	cout << "End find log" << endl << endl;
// 	return NULL;
// }


template<class DATA>
DATA * Bucket<DATA>::find_kmer_linear(kmer_full& kmer, const int64_t begin, const int64_t end) {
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
	// if ((uint)kmer.kmer_s == 2214592512) {
	// if (*(((uint*)&(kmer.kmer_s))+1) == 553648128) {
	// 	debug = true;
	// }
	// if (debug)
		// for (SKL & skmer : skml){
		// 	skmer.print(nucleotides_reserved_memory + params->allocated_bytes * skmer.idx, (kint)34603008UL, *params); cout << endl;
		// }
	if (sorted_size > 0) {
		DATA * ptr = find_kmer_log(kmer);
		// if (debug) {
		// 	print_kmer(kmer.kmer_s, 31); cout << endl;
		// }
		debug = false;
		if (ptr != NULL) {
			return ptr;
		}
	}
	debug = false;

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
