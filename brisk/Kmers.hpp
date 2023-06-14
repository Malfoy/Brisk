#include <emmintrin.h>
#include <tmmintrin.h>
#include <cstdint>
#include <string>
#include <vector>
#include <utility>
#include "parameters.hpp"
#include "hashing.hpp"

using namespace std;



#ifndef KMERS_H
#define KMERS_H


// --- Kint definitions ---

typedef __uint128_t kint;
typedef __uint128_t skint;
//~ typedef uint64_t kint;



// uint64_t hash64shift(uint64_t key);
// static constexpr double precision2bit = std::numeric_limits<uint16_t>::max() * 1.88416938536372; // * log(2)/exp(-1)



// Hash function for kint in robin_hood
// namespace robin_hood {
//  template <> struct hash<kint>
//   {
//     size_t operator()(const kint & x) const
//     {
//       // your code here, e.g. "return hash<int>()(x.value);" 
//       return ((hash64shift(x)) ^ hash64shift(x>>64));
//     }
//   };
// }



// ----- Kmer classes -----
class kmer_full {
public:
	kint kmer_s;
	kint minimizer;
	uint8_t minimizer_idx;
	bool multi_mini;
	vector<int8_t> interleaved;
	DecyclingSet* dede;
	static uint16_t* occ2mer_entropy;

	kmer_full(kint value, uint8_t minimizer_idx, uint8_t minimizer_size, bool multiple_mini,DecyclingSet* dede);
	kmer_full();
	kmer_full(kmer_full&& kmer);
	kmer_full(const kmer_full&& kmer);
	kmer_full & operator=(kmer_full&& kmer);	
	void compute_mini(uint8_t mini_size);
	void print(uint8_t k, uint8_t m) const;
	kint get_compacted(uint8_t m, uint8_t mini_idx)const ;
	// uint64_t get_minimizer() const;
	bool contains_multi_minimizer() const;
	uint8_t prefix_size(const uint8_t k, const uint8_t m) const;
	uint8_t suffix_size() const;
	vector<int> compute_interleaved(const Parameters & params) const;
	// int8_t interleaved_nucleotide(const uint8_t nucl_idx, const uint8_t k, const uint8_t m, bool debug);
	kint get_unhash_kmer_value(uint8_t m) const;
	kmer_full hash_kmer_minimizer_copy(uint8_t m) const;
	void hash_kmer_minimizer_inplace(uint8_t m);
	void unhash_kmer_minimizer(uint8_t m);
	void copy(const kmer_full& kmer);
private:
	void replace_slice(kint replacement, size_t position, size_t length);
};



class SuperKmerEnumerator {
public:
	SuperKmerEnumerator(string & s, const uint8_t k, const uint8_t m,DecyclingSet* dede);
	kint next(vector<kmer_full> & kmers);

private:
	// Sequence and position in it
	string& seq;
	uint64_t seq_idx;
	DecyclingSet* dede;

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
	uint64_t mini_hash;
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
			cerr<<"T";
		}
		if(nuc==3){
			cerr<<"G";
		}
		if(nuc==1){
			cerr<<"C";
		}
		if(nuc==0){
			cerr<<"A";
		}
		if (nuc>=4){
			cerr<<nuc<<endl;
			cerr<<"WTF"<<endl;
		}
		anc>>=2;
	}
	cerr<<endl;
}


string kmer2str(kint num, uint k);
kint str2num(const std::string& str);
kmer_full * str2kmer(const std::string & str, const uint8_t m);
// Return the canonical minimizer for a uint64 sequence.
uint64_t get_minimizer(kint seq, const uint8_t k, uint8_t& min_position, const uint8_t m, bool & reversed, bool & multiple,const uint64_t,DecyclingSet* dede);
string getCanonical(const string& str);
// void string_to_kmers_by_minimizer(string & seq, vector<vector<kmer_full> > & kmers, uint8_t k, uint8_t m);
kint string_to_kmers_by_minimizer(string & seq, vector<kmer_full> & kmers, const uint8_t k, const uint8_t m);



#endif
