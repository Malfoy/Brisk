#include <stdint.h>
#include <string>
#include <assert.h>
#include <vector>

#include <iostream>

#include "Kmers.hpp"
#include "DenseMenuYo.hpp"

using namespace std;


template <class DATA>
class Brisk {
private:
	DenseMenuYo<DATA> * menu;

	uint64_t getMemorySelfMaxUsed() const;
public:
	uint8_t k;
	uint8_t m;

	Brisk(uint8_t k, uint8_t m);

	DATA * insert(kmer_full & kmer);
	DATA * get(kmer_full & kmer) const;

	void protect_data(const kmer_full & kmer);
	void unprotect_data(const kmer_full & kmer);

	bool next(kmer_full & kmer);
	void restart_kmer_enumeration();
	void stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & nb_cursed, uint64_t & memory_usage) const;
};


template<class DATA>
Brisk<DATA>::Brisk(uint8_t k, uint8_t m) {
	this->k = k;
	this->m = m;

	this->menu = new DenseMenuYo<DATA>(k, m);

	assert(m % 2 == 1);
	assert(m < k);
}


template<class DATA>
DATA * Brisk<DATA>::get(kmer_full & kmer) const {
	return this->menu->get_kmer(kmer);
}

template <class DATA>
DATA * Brisk<DATA>::insert(kmer_full & kmer) {
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
void Brisk<DATA>::stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers, uint64_t & nb_cursed, uint64_t & memory_usage) const {
	memory_usage = this->getMemorySelfMaxUsed();
	return this->menu->stats(nb_buckets, nb_skmers, nb_kmers, nb_cursed);
}
