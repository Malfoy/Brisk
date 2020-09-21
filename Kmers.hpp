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
	kint prefix;
	kint suffix;
	bool multi_mini;

	kmer_full(kint value, uint8_t minimizer_idx, uint8_t m, bool multiple_mini);
	void print(uint8_t k, uint8_t m) const;
	kint get_compacted()const ;
	// uint64_t get_minimizer() const;
	bool contains_multi_minimizer() const;
};



// ----- Usefull binary kmer functions -----

void print_kmer(kint num,uint64_t n);
string kmer2str(kint num, uint k);
kint str2num(const std::string& str);
// Return the canonical minimizer for a uint64 sequence.
int64_t get_minimizer(kint seq, uint8_t k, int8_t& position, uint8_t m);
string getCanonical(const string& str);

// void string_to_kmers_by_minimizer(string & seq, vector<vector<kmer_full> > & kmers, uint8_t k, uint8_t m);
kint string_to_kmers_by_minimizer(string & seq, vector<kmer_full> & kmers, const uint8_t k, const uint8_t m);




#endif
