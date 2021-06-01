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


	DATA * insert(kmer_full & kmer);
	DATA * get(kmer_full & kmer);
	
	vector<DATA *> insert_superkmer( vector<kmer_full>& v, vector<bool>& newly_inserted);
	vector<DATA *> get_superkmer( vector<kmer_full>& v);

	vector<DATA *> insert_sequence(const string& str);
	vector<DATA *> get_sequence(const string& str);

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
DATA * Brisk<DATA>::get(kmer_full & kmer) {
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
	return this->menu->get_kmer(kmer);
}


template<class DATA>
vector<DATA *> Brisk<DATA>::get_sequence(const string& str) {
	vector<DATA *> result;
	// Line too short
	if (str.size() < this.params.k){
		return result;
	}
	vector<kmer_full> superkmer;
	SuperKmerEnumerator enumerator(str, this.params.k, this.params.m);
	enumerator.next(superkmer);
	while (superkmer.size() > 0) {
		// Add the values
		for (kmer_full & kmer : superkmer) {
			result.push_back(this->menu->get_kmer(kmer));
		}

	}
	return result;
}



template<class DATA>
vector<DATA *> Brisk<DATA>::get_superkmer( vector<kmer_full>& superkmer) {
	vector<DATA *> result;
	if (superkmer.size() > 0) {
		// for (kmer_full & kmer : superkmer) {
		// 	result.push_back(this->menu->get_kmer(kmer));
		// }
		return this->menu->get_kmer_vector(superkmer);

	}
	return result;
}


template<class DATA>
vector<DATA *> Brisk<DATA>::insert_sequence(const string& str) {
	vector<DATA *> result;
	// Line too short
	if (str.size() < this.params.k){
		return result;
	}
	vector<kmer_full> superkmer;
	SuperKmerEnumerator enumerator(str, this.params.k, this.params.m);
	enumerator.next(superkmer);
	while (superkmer.size() > 0) {
		// Add the values
		uint64_t small_minimizer = (uint32_t)(superkmer[0].minimizer & this->menu->mini_reduc_mask);
		uint32_t mutex_idx = (small_minimizer%this->menu->mutex_number);
		omp_set_lock(this->menu->MutexBucket[mutex_idx]);
		for (kmer_full & kmer : superkmer) {
			result.push_back(this->menu->insert_kmer_no_mutex(kmer));
		}
		omp_unset_lock(this->menu->MutexBucket[mutex_idx]);

	}
	return result;
}


template<class DATA>
vector<DATA *> Brisk<DATA>::insert_superkmer(vector<kmer_full>& superkmer, vector<bool>& newly_inserted){
	vector<DATA *> result;
	if (superkmer.size() > 0) {
		// Add the values
		uint64_t minimizer(superkmer[0].minimizer);
		uint64_t small_minimizer = (uint32_t)(minimizer & this->menu->mini_reduc_mask);
		uint32_t mutex_idx = (small_minimizer%this->menu->mutex_number);
		omp_set_lock(&(this->menu->MutexBucket[mutex_idx]));
		result=this->menu->insert_kmer_vector(superkmer,newly_inserted);
		// for (kmer_full & kmer : superkmer) {
		// 	bool newly_inserted_element;
		// 	DATA* lol=this->menu->insert_kmer_no_mutex(kmer,newly_inserted_element);
		// 	newly_inserted.push_back(newly_inserted_element);
		// }
		// for (kmer_full & kmer : superkmer) {
		// 	DATA* lol=this->menu->get_kmer_no_mutex(kmer);
		// 	result.push_back(lol);
		// }
		result=this->menu->get_kmer_vector(superkmer);
		// cout<<"kmer get"<<endl;
		omp_unset_lock(&(this->menu->MutexBucket[mutex_idx]));
	}
	return result;
}






template <class DATA>
DATA * Brisk<DATA>::insert(kmer_full & kmer) {
	#ifdef TIME_ANALYSIS
	#pragma omp critical
	{
		nb_insert += 1;

		if (nb_get + nb_insert == 1000000) {
			auto current_time = std::chrono::system_clock::now();
			chrono::duration<double> elapsed_seconds = current_time - previous_time;
			cout << "1 000 000 ops: " << elapsed_seconds.count() << endl;
			
			previous_time = current_time;
			nb_get = nb_insert = 0;
		}
	}
	#endif
	return this->menu->insert_kmer(kmer);
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
	return this->menu->next(kmer);
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
