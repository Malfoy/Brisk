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
	bool cleared;

	Bucket(Parameters * params);
	Bucket(const Bucket<DATA> &) = delete;
	Bucket(Bucket<DATA> && bucket);
	~Bucket();

	Bucket<DATA> & operator=(const Bucket<DATA>&) = delete;
	Bucket<DATA> & operator=(Bucket<DATA>&& bucket);

	DATA * insert_kmer(const kmer_full & kmer);
	DATA * find_kmer(const kmer_full& kmer);
	vector<DATA *> find_kmer_vector(const vector<kmer_full>& kmers);

	void init_enumerator();
	void next_kmer(kmer_full & kmer, kint minimizer);
	bool has_next_kmer();

	void clear();
	void print();

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

	DATA * find_kmer_unsorted(const kmer_full& kmer);
	// DATA * find_kmer_from_interleave(kmer_full& kmer, SKL & mockskm, uint8_t * mock_nucleotides);
	DATA * find_kmer_linear(const kmer_full& kmer, const int64_t begin, const int64_t end);
	vector<DATA *> find_kmer_linear_vector(const vector<kmer_full>& kmers,vector<DATA *>& result, const uint64_t begin, const uint64_t end);
	DATA * find_kmer_linear_sorted_stop(const kmer_full& kmer, const int64_t begin, const int64_t end);
	vector<DATA*>  find_kmer_linear_sorted_stop_vector(const vector<kmer_full>& kmers,  vector<uint64_t>& begins, const uint64_t end);
	DATA * find_kmer_log(const kmer_full & kmer);
	DATA * find_kmer_log_simple(const kmer_full & kmer);
	vector<DATA *> find_kmer_log_simple_vector(const vector<kmer_full> & kmers);
	DATA * insert_kmer_buffer(const kmer_full & kmer);
	void discard_last_kmer();
	void insert_buffer();
	void data_space_update();
};



template <class DATA>
Bucket<DATA>::Bucket(Parameters * params) {
	this->sorted_size = 0;
	this->buffered_skmer = (SKL *)NULL;
	this->buffered_get = ((kint)1) << (sizeof(kint) * 8 - 1);
	this->buffered_data = 0;
	this->cleared = false;

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
, cleared( bucket.cleared )
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
void Bucket<DATA>::clear() {
	free(this->nucleotides_reserved_memory);
	this->nucleotides_reserved_memory = nullptr;
	free(this->data_reserved_memory);
	this->data_reserved_memory = nullptr;
	this->cleared = true;
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
		auto factor = this->data_reserved_number > 500 ? 1 : 2;
		size_t to_reserve = (int)(factor * this->data_reserved_number);
		this->data_reserved_memory = (DATA *)realloc(this->data_reserved_memory, sizeof(DATA) * (this->data_reserved_number + to_reserve));
		this->data_reserved_number += to_reserve;
		
		this->buffered_data = NULL;
		this->buffered_get = ((kint)1) << (sizeof(kint) * 8 - 1);
	}
}



template <class DATA>
DATA * Bucket<DATA>::insert_kmer(const kmer_full & kmer) {
	// cout << "bucket insert " << kmer2str(kmer.kmer_s, this->params->k) << endl;
	// 0 - Update space for DATA if needed
	this->data_space_update();
	// 1 - Try to compact with the last kmer
	if (buffered_skmer != NULL) {
		// TODO: small minimizer idx ok ?
		bool is_compacted = buffered_skmer->compact_right(
				kmer,
				this->nucleotides_reserved_memory + this->params->allocated_bytes * buffered_skmer->idx,
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
	
	// 3 - Create a new skmer
	DATA * value = this->insert_kmer_buffer(kmer);
	this->nb_kmers += 1;
	this->next_data += 1;
	this->buffered_data = value;
	this->buffered_get = kmer.kmer_s;

	// 2 - Sort if needed
	//TAMPON
	if(skml.size()-sorted_size>100 and skml.size()-sorted_size>sorted_size*1){
		this->insert_buffer();
	}

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
		// a.print(
		// 	nucleotides_reserved_memory + a.idx * params->allocated_bytes,
		// 	1,
		// 	*params
		// ); cout << endl;
		// b.print(
		// 	nucleotides_reserved_memory + b.idx * params->allocated_bytes,
		// 	1,
		// 	*params
		// ); cout << endl;
		// cout << val << endl;
		return val;
	};
	// TODO: We should not sort the last kmer
	sort(skml.begin()+sorted_size, skml.end(), comp_function);
	inplace_merge(skml.begin(), skml.begin()+sorted_size, skml.end(), comp_function);

	buffered_skmer = NULL;
	this->buffered_data = NULL;
	this->buffered_get = ((kint)1) << (sizeof(kint) * 8 - 1);
	sorted_size=skml.size();
}



template <class DATA>
DATA * Bucket<DATA>::insert_kmer_buffer(const kmer_full & kmer){
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
		kmer.get_compacted(params->m_small,
						   kmer.minimizer_idx + (this->params->m_reduc + 1) / 2),
		(int)kmer.minimizer_idx + (this->params->m_reduc + 1) / 2,
		skml.size(),
		this->nucleotides_reserved_memory + (this->skml.size() * params->allocated_bytes),
		this->next_data,
		*params
	);
	buffered_skmer = &(skml[skml.size()-1]);
	value_pointer = this->data_reserved_memory + this->next_data;
	
	return value_pointer;
}



template <class DATA>
void Bucket<DATA>::discard_last_kmer(){
	skml.pop_back();
	buffered_skmer = &(skml[skml.size()-1]);
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
DATA * Bucket<DATA>::find_kmer_log_simple(const kmer_full & kmer) {
	// return (find_kmer_linear_sorted_stop(kmer,0,sorted_size-1));
	// kmer.minimizer_idx+=params->m-params->m_small;
	insert_kmer_buffer(kmer);
	auto comp_function = [&](const SKL & a, const SKL & b) {
		bool val = a.inf(
				nucleotides_reserved_memory + a.idx * params->allocated_bytes,
				b,
				nucleotides_reserved_memory + b.idx * params->allocated_bytes,
				*params
		);
		return val;
	};

	uint64_t bottom(lower_bound(skml.begin(), skml.begin()+sorted_size,skml[skml.size()-1] ,comp_function)-skml.begin());
	discard_last_kmer();
	DATA* search=(find_kmer_linear_sorted_stop(kmer,bottom,sorted_size-1));
	return search;
}



template<class DATA>
vector<DATA *> Bucket<DATA>::find_kmer_log_simple_vector(const vector<kmer_full> & kmers) {
	auto comp_function = [&](const SKL & a, const SKL & b) {
		bool val = a.inf(
				nucleotides_reserved_memory + a.idx * params->allocated_bytes,
				b,
				nucleotides_reserved_memory + b.idx * params->allocated_bytes,
				*params
		);
		return val;
	};
	
	vector<uint64_t> begins(kmers.size());
	for(uint i(0);i<kmers.size(); ++i){
		insert_kmer_buffer(kmers[i]);
		begins[i]= lower_bound(skml.begin(), skml.begin()+sorted_size,skml[skml.size()-1] ,comp_function)-skml.begin();
		discard_last_kmer();
	}
	return find_kmer_linear_sorted_stop_vector(kmers,begins,sorted_size-1);
}



template<class DATA>
DATA * Bucket<DATA>::find_kmer_log(const kmer_full & kmer) {

	vector<int> kmer_interleved = kmer.compute_interleaved(*params);

	vector<int> current_interleved(kmer_interleved);


	// Prepare the search parameters
	int begin = 0; int end = sorted_size - 1;
	vector<int> heap;
	uint interleaved_size = kmer_interleved.size();
	vector<int> middle_interleved(interleaved_size, -2);
	// int * middle_interleved = new int[interleaved_size / 8 + 1];
	int prev_middle = -1;

	// Log search
	while (begin <= end or heap.size() > 0) {
		// Base case - not found, restore previous search from heap
		if (begin > end) {
			uint interleaved_idx = heap.back(); heap.pop_back();
			end = heap.back(); heap.pop_back();
			begin = heap.back(); heap.pop_back();

			// First time modified
			if (current_interleved[interleaved_idx] == -2) {
				// Set all unknown values to -1
				for (uint idx=interleaved_idx ; idx<interleaved_size ; idx += 2)
					current_interleved[idx] = -1;
			}
			// Restore the unknown following -2 and increment the current interleave
			else if (current_interleved[interleaved_idx] == -1) {
				current_interleved[interleaved_idx] = 0;
				for (uint idx=interleaved_idx+2 ; idx<interleaved_size ; idx += 2)
					current_interleved[idx] = -2;
			}
			// In other cases, change the value of the interleaved
			else {
				current_interleved[interleaved_idx] += 1;
				for (uint idx=interleaved_idx+2 ; idx<interleaved_size ; idx += 2)
					if (kmer_interleved[idx] == -2)
						current_interleved[idx] = -2;
			}

			// If remaining possible values, save the context
			if (current_interleved[interleaved_idx] < 3) {
				heap.push_back(begin); heap.push_back(end); heap.push_back(interleaved_idx);
			}

		}

		// Log case - get the middle superkmer
		int middle = begin + (end - begin) / 2;
		if (middle != prev_middle) {
			SKL & mid_skmer = skml[middle];
			mid_skmer.compute_interleaved(middle_interleved, nucleotides_reserved_memory + params->allocated_bytes * mid_skmer.idx, *params);
			prev_middle = middle;
		}


		// Compare middle superkmer and searched superkmer
		bool found = true;
		for (uint idx=0 ; idx<interleaved_size ; idx++) {
			// Undefined current nucleotide
			if (current_interleved[idx] == -2) {
				// Save context for restoration
				heap.push_back(begin); heap.push_back(end); heap.push_back((int)idx);
				// Nothing to do here
				end = begin - 1;
				found = false;
				break;
			}
			else if (current_interleved[idx] != middle_interleved[idx]) {
				found = false;
				if (current_interleved[idx] < middle_interleved[idx])
					end = middle - 1;
				else
					begin = middle + 1;
				break;
			}

		}

		if (found) {
			uint kmer_position = skml[middle].prefix_size(*params) - kmer.prefix_size(params->k, params->m_small);
			return this->data_reserved_memory + skml[middle].data_idx + kmer_position;
		}
	}

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
			uint64_t mini_reduc_offset = (this->params->m_reduc + 1) / 2;
			uint64_t kmer_mini_idx = kmer.minimizer_idx + mini_reduc_offset;
			uint8_t kmer_position = skml[i].size - (skml[i].minimizer_idx - kmer_mini_idx) - 1;
			buffered_data = this->data_reserved_memory + skml[i].data_idx + kmer_position;
			return this->data_reserved_memory + skml[i].data_idx + kmer_position;
		}
	}
	return NULL;
}



template<class DATA>
vector<DATA *> Bucket<DATA>::find_kmer_linear_vector(const vector<kmer_full>& kmers,vector<DATA *>& result, const uint64_t begin, const uint64_t end) {
	for (uint64_t i=begin ; i<=end ; i++) {
		for(uint64_t j(0);j<kmers.size(); ++j){
			if(result[j]==NULL){
				// cout << "i=" << i << " idx=" << (uint64_t)skml[i].idx << " allocated=" << params->allocated_bytes << endl;
				bool is_present = skml[i].is_kmer_present(
					kmers[j],
					this->nucleotides_reserved_memory + params->allocated_bytes * skml[i].idx,
					*params);
				// cout << is_present << endl;
				// cout << "-----" << endl;
				if (is_present) {
					uint64_t mini_reduc_offset = (this->params->m_reduc + 1) / 2;
					uint64_t kmer_mini_idx = kmers[j].minimizer_idx + mini_reduc_offset;
					uint8_t kmer_position = skml[i].size - (skml[i].minimizer_idx - kmer_mini_idx) - 1;
					buffered_data = this->data_reserved_memory + skml[i].data_idx + kmer_position;
					result[j]=this->data_reserved_memory + skml[i].data_idx + kmer_position;
				}
			}
		}
	}
	return result;
}



template<class DATA>
DATA * Bucket<DATA>::find_kmer_linear_sorted_stop(const kmer_full& kmer, const int64_t begin, const int64_t end) {
	vector<int> kmer_interleave=kmer.compute_interleaved(*params);
	for (int i=begin ; i<=end ; i++) {
		bool inferior,superior,equal;
		if(skml[i].kmer_comparison(kmer,kmer_interleave,this->nucleotides_reserved_memory + params->allocated_bytes * skml[i].idx,*params,superior,inferior,equal)){
			if(equal){
				uint64_t mini_reduc_offset = (this->params->m_reduc + 1) / 2;
				uint64_t kmer_mini_idx = kmer.minimizer_idx + mini_reduc_offset;
				uint8_t kmer_position = skml[i].size - (skml[i].minimizer_idx - kmer_mini_idx) - 1;
				buffered_data = this->data_reserved_memory + skml[i].data_idx + kmer_position;
				return this->data_reserved_memory + skml[i].data_idx + kmer_position;
			}else if(superior){
				return NULL;
			}
		}
	}
	return NULL;
}



template<class DATA>
vector<DATA *> Bucket<DATA>::find_kmer_linear_sorted_stop_vector(const vector<kmer_full>& kmers,  vector<uint64_t>& begins, const uint64_t end) {

	vector<vector<int>> kmer_interleaves(kmers.size());
	vector<int> superkmer_interleave_buffer(2* params->k,-3);
	vector<DATA*> result(kmers.size(),NULL);
	for(uint64_t i(0);i<kmers.size(); ++i){
		kmer_interleaves[i]= kmers[i].compute_interleaved(*params);
	}
	// DEBUG TODO
	uint64_t min_value(*min_element(begins.begin(),begins.end()));
	uint64_t min_kmer_indice(0);
	//FOREACH SUPERKMER in RANGE
	for (uint64_t i=min_value; i<=end ; i=max(i+1,min_value)) {
		superkmer_interleave_buffer.assign(2* params->k,-3);
		//FOREACH KMER
		for(uint64_t j(min_kmer_indice);j<kmers.size(); ++j){
			if(begins[j]<=i){
				bool inferior,superior,equal;
				skml[i].kmer_comparison(kmers[j],superkmer_interleave_buffer,kmer_interleaves[j],this->nucleotides_reserved_memory + params->allocated_bytes * skml[i].idx,*params,superior,inferior,equal);
				if(equal){
					uint64_t mini_reduc_offset = (this->params->m_reduc + 1) / 2;
					uint64_t kmer_mini_idx = kmers[j].minimizer_idx + mini_reduc_offset;
					uint8_t kmer_position = skml[i].size - (skml[i].minimizer_idx - kmer_mini_idx) - 1;
					buffered_data = this->data_reserved_memory + skml[i].data_idx + kmer_position;
					result[j]=this->data_reserved_memory + skml[i].data_idx + kmer_position;
					if(min_kmer_indice==j){
						min_kmer_indice++;
					}
					begins[j]=end+1;
					if(begins[j]==min_value){
						min_value =*min_element(begins.begin(),begins.end());
					}
					begins[j]=end+1;
				}else if(superior){
					if(min_kmer_indice==j){
						min_kmer_indice++;
					}
					begins[j]=end+1;
					if(begins[j]==min_value){
						min_value =*min_element(begins.begin(),begins.end());
					}
					
				}
			}
		}
	}
	return result;
}



template <class DATA>
DATA * Bucket<DATA>::find_kmer_unsorted(const kmer_full& kmer) {
	return find_kmer_linear(kmer, sorted_size, skml.size()-1);
}



template <class DATA>
DATA * Bucket<DATA>::find_kmer(const kmer_full& kmer) {
	if (sorted_size > 0) {
		DATA * ptr = find_kmer_log_simple(kmer);
		if (ptr != NULL) {
			return ptr;
		}else{	
		}
	}
	return find_kmer_unsorted(kmer);
}



template <class DATA>
vector<DATA *> Bucket<DATA>::find_kmer_vector(const vector<kmer_full>& kmers) {
	cout << "--- find_kmer_vector ---" << endl;
	for (SKL & sk: skml) {
		sk.print(nucleotides_reserved_memory + params->allocated_bytes * sk.idx, 0b01111000, *params); cout << endl;
	}
	cout << "--- find_kmer_vector ---" << endl;
	vector<DATA *> result(kmers.size(),NULL);
	if(sorted_size > 0){
		cout << sorted_size << endl;
		result = find_kmer_log_simple_vector(kmers);
	}
	return find_kmer_linear_vector(kmers,result,sorted_size,skml.size()-1);
}



template <class DATA>
void Bucket<DATA>::init_enumerator() {
	this->enumeration_skmer_idx = 0;
	this->enumeration_kmer_idx = 0;
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
	if (not has_next_kmer()){
		return;
	}

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
