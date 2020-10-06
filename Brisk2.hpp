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
public:
	uint8_t k;
	uint8_t m;

	Brisk(uint8_t k, uint8_t m);

	DATA * insert(kmer_full & kmer, const kint minimizer);
	DATA * get(kmer_full & kmer, const kint minimizer);
	bool next(kmer_full & kmer);
	void restart_kmer_enumeration();
	void stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers) const;
};


template<class DATA>
Brisk<DATA>::Brisk(uint8_t k, uint8_t m){
	this->k = k;
	this->m = m;

	this->menu = new DenseMenuYo<DATA>(k, m);

	assert(m % 2 == 1);
	assert(m < k);
}


template<class DATA>
DATA * Brisk<DATA>::get(kmer_full & kmer, const kint minimizer) {
	return this->menu->get_kmer(kmer, minimizer);
}

template <class DATA>
DATA * Brisk<DATA>::insert(kmer_full & kmer, const kint minimizer) {
	return this->menu->insert_kmer(kmer, minimizer);
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
void Brisk<DATA>::stats(uint64_t & nb_buckets, uint64_t & nb_skmers, uint64_t & nb_kmers) const {
	return this->menu->stats(nb_buckets, nb_skmers, nb_kmers);
}
