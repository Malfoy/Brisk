#include <emmintrin.h>
#include <tmmintrin.h>
#include <cstdint>
#include <string>
#include "robin_hood.h"




using namespace std;



#ifndef KMERS_H
#define KMERS_H


typedef __uint128_t kint;
typedef __uint128_t skint;
//~ typedef uint64_t kint;



extern robin_hood::unordered_flat_map<string, uint64_t> real_count;
extern uint64_t counting_errors;
extern const uint64_t k;
extern const kint k_mask;
extern const kint compact_mask;
extern const  uint64_t minimizer_size;
extern const  uint64_t super_minimizer_size;
extern const uint64_t min_mask;
extern const uint64_t compacted_size;
extern const uint64_t byte_nuc;






// ----- Usefull binary kmer functions -----

string intToString(kint n);

void print_kmer(kint num,uint64_t n);
string kmer2str(kint num, uint k);
kint str2num(const std::string& str);
// RC functions
uint64_t rcb(const uint64_t& in);
__uint128_t rcb(const __uint128_t& in);
uint64_t rcbc(uint64_t in, uint64_t n);
// Hash functions (TODO: move them)
__m128i mm_bitshift_right(__m128i x, unsigned count);
__m128i mm_bitshift_left(__m128i x, unsigned count);
uint64_t hash64shift(uint64_t key);
// Return the canonical minimizer for a uint64 sequence.
int64_t get_minimizer(kint seq, int8_t& position);
uint64_t reversebits(uint64_t b);
string getCanonical(const string& str);




// ----- Kmer class -----
class kmer_full {
public:
	int8_t minimizer_idx;
	kint kmer_s;
	kint prefix;
	kint suffix;
	kmer_full(int8_t minimizer_idx, kint value);
	kint get_compacted()const ;
	/** Return the minimizer regarding the minimizer_idx property
		* Warning: The minimizer should be canon
		*/
	uint8_t get_minimizer_idx() const;
	// uint64_t get_minimizer() const;
	bool contains_multi_minimizer() const;
};

#endif
