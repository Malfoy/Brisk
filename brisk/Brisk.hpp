#include <stdint.h>
#include <string>
#include <assert.h>
#include <vector>
#include <iostream>
#include <cstdint>
#include "Kmers.hpp"
#include "buckets.hpp"
#include "DenseMenuYo.hpp"
#include "parameters.hpp"



using namespace std;



#ifndef BRISK_H
#define BRISK_H



template <class DATA>
class Brisk {
public:
	uint64_t getMemorySelfMaxUsed() const;
	vector<DATA *> insert_sequence(const string& str, vector<bool>& newly_inserted);
	vector<DATA *> get_sequence(const string& str);	
	DenseMenuYo<DATA> * menu;
	Parameters params;
	Brisk(Parameters & parameters);
	~Brisk();
	DATA * get(kmer_full & kmer) const;
	vector<DATA *> insert_superkmer(vector<kmer_full>& v, vector<bool>& newly_inserted);
	vector<DATA *> get_superkmer(vector<kmer_full>& v);
	void protect_data(const kmer_full & kmer);
	void unprotect_data(const kmer_full & kmer);
	bool next(kmer_full & kmer);
	void restart_kmer_enumeration();
	void stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & memory_usage, uint64_t & largest_bucket) const;
	void reallocate();
};



template<class DATA>
Brisk<DATA>::Brisk(Parameters & parameters): params( parameters ) {
	Bucket<DATA>::setParameters(parameters);
	this->menu = new DenseMenuYo<DATA>(params);
	assert(params.m % 2 == 1);
	assert(params.m < params.k);
}



template<class DATA>
Brisk<DATA>::~Brisk() {
	delete this->menu;
}



template<class DATA>
DATA * Brisk<DATA>::get(kmer_full & kmer) const {
	kmer.hash_kmer_minimizer_inplace(this->params.m);
	DATA * val = this->menu->get_kmer(kmer);
	kmer.unhash_kmer_minimizer(this->params.m);
	return val;
}



void hash_skmer(vector<kmer_full> & skmer, size_t m) {
	for (kmer_full & kmer : skmer) {
		kmer.hash_kmer_minimizer_inplace(m);
	}
}



void unhash_skmer(vector<kmer_full> & skmer, size_t m) {
	for (kmer_full & kmer : skmer)
		kmer.unhash_kmer_minimizer(m);
}



kmer_full update_kmer(const kmer_full& old_kmer, uint8_t k, uint8_t m, DecyclingSet* dede) {
	kint km_val = old_kmer.kmer_s;
	uint8_t min_pos;
	bool reversed;
	get_minimizer(km_val, k, min_pos, m, reversed, ((kint)1<<(2*m))-1,dede);
	if (not reversed)
		return kmer_full(km_val, min_pos, m, dede);
	else
		return kmer_full(km_val, k - m - min_pos, m, dede);
}



template<class DATA>
vector<DATA *> Brisk<DATA>::get_superkmer( vector<kmer_full>& superkmer) {
	vector<DATA *> result;
	if (superkmer.size() > 0) {
		hash_skmer(superkmer, this->params.m);
		// Remove the minimizer suffix
		uint64_t small_minimizer = superkmer[0].minimizer >> (2 * ((this->params.m_reduc + 1) / 2));
		// Remove the minimizer prefix
		uint32_t mutex_idx = ((small_minimizer)%this->menu->mutex_number);

		small_minimizer &= (((kint)1) << (this->params.b * 2)) - 1;
		omp_set_lock(&(this->menu->MutexBucket[mutex_idx]));
		result=this->menu->get_kmer_vector(superkmer);
		omp_unset_lock(&(this->menu->MutexBucket[mutex_idx]));
		unhash_skmer(superkmer, this->params.m);
	}
	return result;
}



template<class DATA>
vector<DATA *> Brisk<DATA>::insert_superkmer(vector<kmer_full>& superkmer, vector<bool>& newly_inserted){
	/*cout << "largest bucket: " << this->menu->largest_bucket << endl;
	if (this->menu->largest_bucket >= 100){
		cout << "doublaj" << endl;
		reallocate();
		cout << "doublage" << endl;
	}*/

	vector<DATA *> result;
	if (superkmer.size() > 0) {
		hash_skmer(superkmer, this->params.m);
		// Remove the minimizer suffix
		uint64_t small_minimizer = superkmer[0].minimizer >> (2 * ((this->params.m_reduc + 1) / 2));
		// Remove the minimizer prefix
		small_minimizer &= (((kint)1) << (this->params.b * 2)) - 1;
		// Add the values
		uint32_t mutex_idx = ((small_minimizer)%this->menu->mutex_number);
		omp_set_lock(&(this->menu->MutexBucket[mutex_idx]));
		this->menu->insert_kmer_vector(superkmer,newly_inserted);
		result=this->menu->get_kmer_vector(superkmer);
		omp_unset_lock(&(this->menu->MutexBucket[mutex_idx]));
		unhash_skmer(superkmer, this->params.m);
	}
	return result;
}



template <class DATA>
void Brisk<DATA>::protect_data(const kmer_full & kmer) {
	this->menu->protect_data(kmer);
}



template <class DATA>
void Brisk<DATA>::unprotect_data(const kmer_full & kmer) {
	this->menu->unprotect_data(kmer);
}



template<class DATA>
bool Brisk<DATA>::next(kmer_full & kmer) {
	bool has_next = this->menu->next(kmer);
	if (has_next) {
		kmer.unhash_kmer_minimizer(this->params.m);
	}
	return has_next;
}



template<class DATA>
void Brisk<DATA>::restart_kmer_enumeration() {
	this->menu->restart_kmer_enumeration();
}



template<class DATA>
uint64_t Brisk<DATA>::getMemorySelfMaxUsed () const{
	uint64_t result = 0;
	struct rusage usage;
	if (getrusage(RUSAGE_SELF, &usage)==0)  {  result = usage.ru_maxrss;  }
	return result;
}



template<class DATA>
void Brisk<DATA>::stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & memory_usage, uint64_t & largest_bucket) const {
	memory_usage = this->getMemorySelfMaxUsed();
	return this->menu->stats(nb_buckets, nb_skmers, nb_kmers, largest_bucket);
}



template<class DATA>
void Brisk<DATA>::reallocate(){
	bool new_insert = false;
	uint32_t bucket_size = 0;
	uint8_t new_m = this->params.m+2;
	uint8_t new_b = this->params.b+2;
	Parameters* new_params = new Parameters(this->params.k,new_m,new_b);
	DenseMenuYo<DATA>* big_brother = new DenseMenuYo<DATA>(*new_params);
	kmer_full old_kmer(0,0, this->params.m, this->params.dede);
	while (this->next(old_kmer)) {
		old_kmer.hash_kmer_minimizer_inplace(this->params.m);
		DATA* old_value = this->menu->get_kmer_no_mutex(old_kmer);
		old_kmer.unhash_kmer_minimizer(this->params.m);	

		kmer_full new_kmer=update_kmer(old_kmer, this->params.k, new_m, new_params->dede);
		new_kmer.hash_kmer_minimizer_inplace(new_m);
		DATA* value = big_brother->insert_kmer_no_mutex(new_kmer, new_insert,bucket_size);
		new_kmer.unhash_kmer_minimizer(new_m);
		*value = *old_value;
	}
	delete this->menu;
	this->menu=big_brother;
	this->params=*new_params;
}



#endif