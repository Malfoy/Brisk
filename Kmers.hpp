#include <emmintrin.h>
#include <tmmintrin.h>
#include <cstdint>
#include <string>
#include <vector>
#include "robin_hood.h"




using namespace std;



#ifndef KMERS_H
#define KMERS_H


// --- Kint definitions ---

typedef __uint128_t kint;
typedef __uint128_t skint;
//~ typedef uint64_t kint;

// Hash function for kint in robin_hood
namespace std {
 template <> struct hash<kint>
  {
    size_t operator()(const kint & x) const
    {
      /* your code here, e.g. "return hash<int>()(x.value);" */
      return ((uint64_t)x+(uint64_t)(x>>64));
    }
  };
}


// ----- Kmer class -----
class kmer_full {
public:
	int8_t minimizer_idx;
	kint kmer_s;
	// kint prefix;
	// kint suffix;
	bool multi_mini;

	kmer_full(kint value, uint8_t minimizer_idx, bool multiple_mini);
	void print(uint8_t k, uint8_t m) const;
	kint get_compacted(uint8_t m)const ;
	// uint64_t get_minimizer() const;
	bool contains_multi_minimizer() const;
};



// ----- Usefull binary kmer functions -----

template<typename T>
void print_kmer(T num, uint8_t n){
	num &= ((T)1 << (2*n)) - 1;
	T anc = (T)1<<(2 * (n - 1));
	for(uint64_t i(0);i<n and anc!=0;++i){
		uint64_t nuc = num/anc;
		num = num % anc;
		if(nuc==2){
			cout<<"T";
		}
		if(nuc==3){
			cout<<"G";
		}
		if(nuc==1){
			cout<<"C";
		}
		if(nuc==0){
			cout<<"A";
		}
		if (nuc>=4){
			cout<<nuc<<endl;
			cout<<"WTF"<<endl;
		}
		anc>>=2;
	}
}
string kmer2str(kint num, uint k);
kint str2num(const std::string& str);
// Return the canonical minimizer for a uint64 sequence.
int64_t get_minimizer(kint seq, uint8_t k, int8_t& position, uint8_t m);
string getCanonical(const string& str);

// void string_to_kmers_by_minimizer(string & seq, vector<vector<kmer_full> > & kmers, uint8_t k, uint8_t m);
kint string_to_kmers_by_minimizer(string & seq, vector<kmer_full> & kmers, const uint8_t k, const uint8_t m);




#endif
