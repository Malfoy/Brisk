#include <emmintrin.h>
#include <tmmintrin.h>
#include <cstdint>
#include <string>
#include <vector>
#include <utility>
#include "robin_hood.h"
#include "sparse_growth_policy.h"
#include "sparse_hash.h"
#include "sparse_map.h"

using namespace std;



#ifndef KMERS_H
#define KMERS_H


// --- Kint definitions ---

typedef __uint128_t kint;
typedef __uint128_t skint;
//~ typedef uint64_t kint;


uint64_t hash64shift(uint64_t key);
uint64_t bfc_hash_64(uint64_t key, uint64_t mask);
uint64_t bfc_hash_64_inv(uint64_t key, uint64_t mask);
static constexpr double precision2bit = std::numeric_limits<uint16_t>::max() * 1.88416938536372; // * log(2)/exp(-1)



// Hash function for kint in robin_hood
namespace robin_hood {
 template <> struct hash<kint>
  {
    size_t operator()(const kint & x) const
    {
      // your code here, e.g. "return hash<int>()(x.value);" 
      return ((hash64shift(x)) ^ hash64shift(x>>64));
    }
  };
}



// ----- Kmer classes -----
class kmer_full {
public:
	kint kmer_s;
	kint minimizer;
	uint8_t minimizer_idx;
	bool multi_mini;
	vector<int8_t> interleaved;
	static uint16_t* occ2mer_entropy;

	kmer_full(kint value, uint8_t minimizer_idx, uint8_t minimizer_size, bool multiple_mini);
	kmer_full();
	kmer_full(kmer_full&& kmer);
	kmer_full(const kmer_full&& kmer);
	kmer_full & operator=(kmer_full&& kmer);	
	void compute_mini(uint8_t mini_size);
	void print(uint8_t k, uint8_t m) const;
	kint get_compacted(uint8_t m)const ;
	// uint64_t get_minimizer() const;
	bool contains_multi_minimizer() const;
	uint8_t prefix_size(const uint8_t k, const uint8_t m) const;
	uint8_t suffix_size() const;
	vector<int> compute_interleaved(const uint8_t k, const uint8_t m) const;
	// int8_t interleaved_nucleotide(const uint8_t nucl_idx, const uint8_t k, const uint8_t m, bool debug);
	void hash_kmer_body(uint8_t m, uint64_t mask_large_minimizer);
	void unhash_kmer_body(uint8_t m, uint64_t mask_large_minimizer);
	kint get_unhash_kmer_body(uint8_t m, uint64_t mask_large_minimizer)const;
	double bimer_entropy(int k );
	void initocc2mer_entropy(int k);
	void copy(const kmer_full& kmer);
};



class SuperKmerEnumerator {
public:
	SuperKmerEnumerator(string & s, const uint8_t k, const uint8_t m);
	kint next(vector<kmer_full> & kmers,bool hash=true);

private:
	// Sequence and position in it
	string& seq;
	uint64_t seq_idx;

	// kmer size and minimizer size
	uint8_t k;
	kint k_mask;
	uint8_t m;
	kint m_mask;

	// Previous kmer read
	bool saved;
	kmer_full saved_kmer;

	// Current values for kmers
	kint current_kmer;
	kint current_rc_kmer;

	// Needed variables
	kint mini_candidate;
	kint rc_mini_candidate;
	bool reversed;
	bool multiple;
	uint8_t mini_pos;
	uint64_t mini;
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


 inline uint64_t bfc_hash_64(uint64_t key, uint64_t mask)
{
	
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}




 inline uint64_t bfc_hash_64_inv(uint64_t key, uint64_t mask)
{
	uint64_t tmp;
 
	// Invert key = key + (key << 31)
	tmp = (key - (key << 31));
	key = (key - (tmp << 31)) & mask;
 
	// Invert key = key ^ (key >> 28)
	tmp = key ^ key >> 28;
	key = key ^ tmp >> 28;
 
	// Invert key *= 21
	key = (key * 14933078535860113213ull) & mask;
 
	// Invert key = key ^ (key >> 14)
	tmp = key ^ key >> 14;
	tmp = key ^ tmp >> 14;
	tmp = key ^ tmp >> 14;
	key = key ^ tmp >> 14;
 
	// Invert key *= 265
	key = (key * 15244667743933553977ull) & mask;
 
	// Invert key = key ^ (key >> 24)
	tmp = key ^ key >> 24;
	key = key ^ tmp >> 24;
 
	// Invert key = (~key) + (key << 21)
	tmp = ~key;
	tmp = ~(key - (tmp << 21));
	tmp = ~(key - (tmp << 21));
	key = ~(key - (tmp << 21)) & mask;
 
	return key;
}


string kmer2str(kint num, uint k);
kint str2num(const std::string& str);
kmer_full * str2kmer(const std::string & str, const uint8_t m);
// Return the canonical minimizer for a uint64 sequence.
uint64_t get_minimizer(kint seq, const uint8_t k, uint8_t& min_position, const uint8_t m, bool & reversed, bool & multiple,const uint64_t);
string getCanonical(const string& str);
// void string_to_kmers_by_minimizer(string & seq, vector<vector<kmer_full> > & kmers, uint8_t k, uint8_t m);
kint string_to_kmers_by_minimizer(string & seq, vector<kmer_full> & kmers, const uint8_t k, const uint8_t m);




#endif
