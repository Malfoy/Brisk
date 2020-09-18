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
	DenseMenuYo * menu;
public:
	uint8_t k;
	uint8_t m;

	Brisk(uint8_t k, uint8_t m);

	DATA * insert(vector<kmer_full> & kmer);
	// DATA * get_pointer(kint kmer);
};


template<class DATA>
Brisk<DATA>::Brisk(uint8_t k, uint8_t m){
	this->k = k;
	this->m = m;

	this->menu = new DenseMenuYo(k, m);

	assert(m % 2 == 1);
	assert(m < k);
}


template <class DATA>
DATA * Brisk<DATA>::insert(vector<kmer_full> & kmer) {
	cout << "TODO: insert" << endl;
	return NULL;
}
