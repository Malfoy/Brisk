#include <stdint.h>
#include <string>
#include <assert.h>
#include <vector>
#include <iostream>
// #define TIME_ANALYSIS
#include "Kmers.hpp"
#include "DenseMenuYo.hpp"
#include "parameters.hpp"



using namespace std;



#ifndef BRISK_H
#define BRISK_H

#ifdef TIME_ANALYSIS
#include <chrono>
#endif



template <class DATA>
class Brisk {
private:

	uint64_t getMemorySelfMaxUsed() const;
	// TODO: Debug before pushing them public. See insert in the code for more details.
	vector<DATA *> insert_sequence(const string& str, vector<bool>& newly_inserted);
	vector<DATA *> get_sequence(const string& str);
public:
	DenseMenuYo<DATA> * menu;
	
#ifdef TIME_ANALYSIS
	uint64_t nb_get;
	uint64_t nb_insert;
	std::chrono::system_clock::time_point previous_time;
#endif
	Parameters params;

	Brisk(Parameters & parameters);
	~Brisk();


	// DATA * insert(kmer_full & kmer);
	DATA * get(kmer_full & kmer) const;
	
	vector<DATA *> insert_superkmer(vector<kmer_full>& v, vector<bool>& newly_inserted);
	vector<DATA *> get_superkmer(vector<kmer_full>& v);

	// vector<DATA *> insert_sequence(const string& str, vector<bool>& newly_inserted);
	// vector<DATA *> get_sequence(const string& str);

	void protect_data(const kmer_full & kmer);
	void unprotect_data(const kmer_full & kmer);

	bool next(kmer_full & kmer);
	void restart_kmer_enumeration();
	void stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & nb_cursed, uint64_t & memory_usage, uint64_t & largest_bucket) const;
};



template<class DATA>
Brisk<DATA>::Brisk(Parameters & parameters): params( parameters ) {
	this->menu = new DenseMenuYo<DATA>(params);

	#ifdef TIME_ANALYSIS
	previous_time = std::chrono::system_clock::now();
	nb_get = nb_insert = 0;
	#endif

	assert(params.m % 2 == 1);
	assert(params.m < params.k);
}



template<class DATA>
Brisk<DATA>::~Brisk() {
	delete this->menu;
}


template<class DATA>
DATA * Brisk<DATA>::get(kmer_full & kmer) const {
	#ifdef TIME_ANALYSIS
	#pragma omp critical
	{
		nb_get += 1;

		if (nb_get + nb_insert == 1000000) {
			auto current_time = std::chrono::system_clock::now();
			chrono::duration<double> elapsed_seconds = current_time - previous_time;
			cout << "1 000 000 ops: " << elapsed_seconds.count() << endl;
			
			previous_time = current_time;
			nb_get = nb_insert = 0;
		}
	}
	#endif
	kmer.hash_kmer_minimizer_inplace(this->params.m);
	DATA * val = this->menu->get_kmer(kmer);
	kmer.unhash_kmer_minimizer(this->params.m);
	return val;
}


void hash_skmer(vector<kmer_full> & skmer, size_t m) {
	// Replace all the minimizers in the kmers
	for (kmer_full & kmer : skmer) {
		kmer.hash_kmer_minimizer_inplace(m);
	}
}


void unhash_skmer(vector<kmer_full> & skmer, size_t m) {
	// Replace all the minimizers in the kmers
	for (kmer_full & kmer : skmer)
		kmer.unhash_kmer_minimizer(m);
}


template<class DATA>
vector<DATA *> Brisk<DATA>::get_superkmer( vector<kmer_full>& superkmer) {
	vector<DATA *> result;
	if (superkmer.size() > 0) {
		hash_skmer(superkmer, this->params.m);
		result = this->menu->get_kmer_vector(superkmer);
		unhash_skmer(superkmer, this->params.m);
	}
	return result;
}



template<class DATA>
vector<DATA *> Brisk<DATA>::insert_superkmer(vector<kmer_full>& superkmer, vector<bool>& newly_inserted){

	vector<DATA *> result;
	if (superkmer.size() > 0) {
		hash_skmer(superkmer, this->params.m);

		// Remove the minimizer suffix
		uint64_t small_minimizer = superkmer[0].minimizer >> (2 * ((this->params.m_reduc + 1) / 2));
		// Remove the minimizer prefix
		small_minimizer &= (((kint)1) << (this->params.m_small * 2)) - 1;

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
void Brisk<DATA>::stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & nb_cursed, uint64_t & memory_usage, uint64_t & largest_bucket) const {
	memory_usage = this->getMemorySelfMaxUsed();
	return this->menu->stats(nb_buckets, nb_skmers, nb_kmers, nb_cursed, largest_bucket);
}


#endif
