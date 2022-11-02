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
	DATA * get(kmer_full & kmer);
	
	vector<DATA *> insert_superkmer(const  vector<kmer_full>& v, vector<bool>& newly_inserted);
	vector<DATA *> get_superkmer(const vector<kmer_full>& v);


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


// TODO: hash kmer minimizer
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
		auto vec=get_superkmer(superkmer);
		result.insert(result.end(),vec.begin(),vec.end());
	}
	return result;
}



template<class DATA>
vector<DATA *> Brisk<DATA>::get_superkmer(const vector<kmer_full>& superkmer) {
	vector<DATA *> result;
	if (superkmer.size() > 0) {
		// Replace all the minimizers in the kmers by their hash values
		vector<kmer_full> hashed_skmer;
		for (const kmer_full & kmer : superkmer) {
			hashed_skmer.emplace_back(kmer.hash_kmer_minimizer_copy(this->params.m));
		}

		return this->menu->get_kmer_vector(hashed_skmer);
	}
	return result;
}



template<class DATA>
// TODO: This function has a bug if multiple skmers share the same minimizer and if there is a reallocation of a bucket for the second add. In this case, the DATA * from the first skmer are not relevant anymore.
vector<DATA *> Brisk<DATA>::insert_sequence(const string& str,vector<bool>& newly_inserted) {
	vector<DATA *> result;
	// Line too short
	if (str.size() < this.params.k){
		return result;
	}
	vector<kmer_full> superkmer;
	SuperKmerEnumerator enumerator(str, this.params.k, this.params.m);
	enumerator.next(superkmer);
	while (superkmer.size() > 0){
		// Add the values
		vector<bool> newly_inserted_local;
		vector<DATA*> vec(this->insert_superkmer(superkmer,newly_inserted));
		superkmer.clear();
		result.insert(result.end(),vec.begin(),vec.end());
		newly_inserted.insert(newly_inserted.end(),newly_inserted_local.begin(),newly_inserted_local.end());
	}//TODO MISSING A NEXT HERE
	return result;
}



template<class DATA>
vector<DATA *> Brisk<DATA>::insert_superkmer(const vector<kmer_full>& superkmer, vector<bool>& newly_inserted){

	vector<DATA *> result;
	if (superkmer.size() > 0) {
		// Replace all the minimizers in the kmers by their hash values
		vector<kmer_full> hashed_skmer;
		for (const kmer_full & kmer : superkmer) {
			hashed_skmer.emplace_back(kmer.hash_kmer_minimizer_copy(this->params.m));
		}

		// TODO: seems like a shifting bug (odd numbers hell)
        uint64_t small_minimizer =  (((hashed_skmer[0].minimizer& this->menu->mini_reduc_mask))>>(params.m-params.m_small));

		// Add the values
		uint32_t mutex_idx = (((small_minimizer))%this->menu->mutex_number);
		omp_set_lock(&(this->menu->MutexBucket[mutex_idx]));
		this->menu->insert_kmer_vector(hashed_skmer,newly_inserted);
		result=this->menu->get_kmer_vector(hashed_skmer);
		omp_unset_lock(&(this->menu->MutexBucket[mutex_idx]));
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
