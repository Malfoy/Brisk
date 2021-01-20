#include <iostream>
#include <algorithm>
#include "Kmers.hpp"
#include "pow2.hpp"
#include <cmath>




using namespace std;



// RC functions
uint64_t rcbc(uint64_t in, uint64_t n);
// Hash functions (TODO: move them)
uint64_t hash64shift(uint64_t key);



// ----- Kmer class -----
kmer_full::kmer_full(kint value, uint8_t minimizer_idx, uint8_t minimizer_size, bool multiple_mini) {
	this->minimizer_idx = minimizer_idx;
	this->kmer_s = value;
	// Shift the kmer to align minizer on the right
	this->minimizer = value >> (2 * minimizer_idx);
	// Mask bits that are not part of the minimizer
	this->minimizer &= ((kint)1 << (2 * minimizer_size))- 1;
	this->multi_mini = multiple_mini;
}



kmer_full::kmer_full(){}



kmer_full::kmer_full(kmer_full&& kmer) :
	kmer_s(kmer.kmer_s),
	minimizer(kmer.minimizer),
	minimizer_idx(kmer.minimizer_idx),
	multi_mini(kmer.multi_mini),
	interleaved(move(interleaved))
{
	// cout << "move constructed" << endl;
}



kmer_full & kmer_full::operator=(kmer_full&& kmer) {
	this->minimizer_idx = kmer.minimizer_idx;
	this->kmer_s = kmer.kmer_s;
	this->minimizer = kmer.minimizer;
	this->multi_mini = kmer.multi_mini;
	this->interleaved = move(kmer.interleaved);
	// cout << "move assigned" << endl;
	return *this;
}



void kmer_full::compute_mini(uint8_t mini_size) {
	this->minimizer = this->kmer_s >> (2 * this->minimizer_idx);
	this->minimizer &= ((kint)1 << (2 * mini_size))- 1;
}



void kmer_full::print(uint8_t k, uint8_t m) const {
	print_kmer(this->kmer_s, k); cout << endl;
	for (uint8_t i=0 ; i<k - m - this->minimizer_idx ; i++)
		cout << " ";
	kint mini = this->kmer_s >> (2 * this->minimizer_idx);
	print_kmer(mini, m); cout << endl;
}



static int const lookup[4][256] = {
	{0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3},
	{0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}
};



vector<int> kmer_full::compute_interleaved(const uint8_t k, const uint8_t m) const {
	vector<int> interleaved;

	uint8_t suff_size = this->suffix_size();
	uint8_t pref_size = this->prefix_size(k, m);
	uint16_t max_idx = 2 * (uint16_t)max(suff_size, pref_size);
	interleaved.resize(max_idx, -2);

	// Compute suffix interleved
	for (uint offset=0 ; offset<suff_size ; offset++) {
		uint8_t nucl_idx = suff_size - 1 - offset;
		uint8_t byte = ((uint8_t*)&kmer_s)[nucl_idx/4];
		interleaved[2*offset] = lookup[nucl_idx%4][byte];
	}

	// Compute prefix interleved
	for (uint offset=0 ; offset<pref_size ; offset++) {
		uint8_t nucl_idx = suff_size + m + offset;
		uint8_t byte = ((uint8_t*)&kmer_s)[nucl_idx/4];
		interleaved[2*offset+1] = lookup[nucl_idx%4][byte];
	}

	return interleaved;
}



uint8_t kmer_full::suffix_size() const {
	return minimizer_idx;
}



uint8_t kmer_full::prefix_size(const uint8_t k, const uint8_t m) const {
	return k - m - minimizer_idx;
}


// SUFFIX IS AT RIGHT!!!!!!DO NOT CHANGE THIS
kint kmer_full::get_compacted(uint8_t m) const {
	kint mask = (((kint)1) << (minimizer_idx * 2)) - 1;
	kint suffix = kmer_s & mask;

	mask = ~mask;
	kint prefix = (kmer_s >> (2 * m));
	prefix = prefix & mask;
	
	return prefix + suffix;
}



bool kmer_full::contains_multi_minimizer() const {
	return this->multi_mini;
}


void kmer_full::hash_kmer_body(uint8_t m, uint32_t mask_large_minimizer){
	//return;
	kint mask = (((kint)1) << (minimizer_idx * 2)) - 1;
	kint suffix = kmer_s & mask;
	kmer_s >>=(2*minimizer_idx);
	mask=(((kint)1) << (m * 2 )) - 1;
	kint new_minimizer=kmer_s & mask;
	kmer_s >>=(2*m);
	kint prefix =kmer_s;
	new_minimizer=bfc_hash_64(minimizer,mask_large_minimizer);
	kmer_s=suffix+(new_minimizer<<(2*minimizer_idx))+(prefix<<(2*(minimizer_idx+m)));
	minimizer=new_minimizer;
}



double kmer_full::bimer_entropy(int k){
	//cout<<"bimer_entropy"<<endl;
	double entropy(0);
	kint kmer=kmer_s;
	uint8_t counts[16] = {0};
	for(int i = 0; i < k-1; i++) {
		counts[kmer & kint(15u)]++;
		kmer >>= 2;
    }
	for(int i(0); i <16;++i){
		entropy+=occ2mer_entropy[counts[i]];
		//cout<<occ2mer_entropy[counts[i]]<<" "<<entropy/precision2bit<<endl;
	}
	//cin.get();
	return entropy/precision2bit;
}

uint16_t* kmer_full::occ2mer_entropy;

void kmer_full::initocc2mer_entropy(int k){
	//cout<<"init"<<endl;
	occ2mer_entropy=new uint16_t[k];
	occ2mer_entropy[0]=0;
	occ2mer_entropy[k-1]=0;
	for(int i(1);i<k-1;++i){
		double p = i / double(k);
		occ2mer_entropy[i]=uint16_t(-log2(p)*p*precision2bit);
		//cout<<i<<" "<<occ2mer_entropy[i]<<endl;
	}
	//cin.get();
}



void kmer_full::unhash_kmer_body(uint8_t m, uint32_t mask_large_minimizer){
	kint mask = (((kint)1) << (minimizer_idx * 2)) - 1;
	kint suffix = kmer_s & mask;
	kmer_s >>=(2*minimizer_idx);
	mask=(((kint)1) << (m * 2 )) - 1;
	kint minimizer=kmer_s & mask;
	kmer_s >>=(2*m);
	kint prefix =kmer_s;
	minimizer=bfc_hash_64_inv(minimizer,mask_large_minimizer);
	kmer_s=suffix+(minimizer<<(2*minimizer_idx))+(prefix<<(2*(minimizer_idx+m)));
}



// ----- Useful binary kmer functions -----



string kmer2str(__uint128_t num, uint k) {
	string res;
	Pow2<__uint128_t> anc(2 * (k - 1));
	for (uint64_t i(0); i < k; ++i) {
		uint64_t nuc = num / anc;
		num             = num % anc;
		if (nuc == 3) {
			res += "G";
		}
		if (nuc == 2) {
			res += "T";
		}
		if (nuc == 1) {
			res += "C";
		}
		if (nuc == 0) {
			res += "A";
		}
		if (nuc >= 4) {
			cout << "WTF kmer2str" << endl;
			cout<<kmer2str(num,k)<<endl;
			// cout<<(uint6)anc.value()<<endl;
			cout<<nuc<<endl;
			return "";
		}
		anc >>= 2;
	}
	return res;
}


kint str2num(const string& str) {
	kint res(0);
	for (uint64_t i(0); i < str.size(); i++) {
		res <<= 2;
		res += (str[i] / 2) % 4;
	}
	return res;
}



kmer_full * str2kmer(const std::string & str, const uint8_t m) {
	kint km_val(str2num(str));

	uint8_t min_pos;
	bool reversed, multiple;
	
	get_minimizer(km_val, str.size(), min_pos, m, reversed, multiple);

	kmer_full * kmer;

	if (not reversed)
		kmer = new kmer_full(km_val, min_pos, m, multiple);
	else
		kmer = new kmer_full(km_val, str.size() - m - min_pos, m, multiple);

	return kmer;
}



kint str2num2(const string& str) {
	kint res(0);
	for (uint64_t i(1); i <= str.size(); i++) {
		res += ((str[i-1] / 2) % 4)<<(str.size()-i);
	}
	return res;
}


uint64_t rcbc(uint64_t in, uint64_t n) {
	assume(n <= 32, "n=%u > 32", n);
	// Complement, swap byte order
	uint64_t res = __builtin_bswap64(in ^ 0xaaaaaaaaaaaaaaaa);
	// Swap nuc order in bytes
	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
	const uint64_t c2 = 0x3333333333333333;
	res               = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
	res               = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc
	// Realign to the right
	res >>= 64 - 2 * n;
	return res;
}



uint64_t canonize(uint64_t x, uint64_t n) {
	return min(x, rcbc(x, n));
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



uint64_t hash64shift(uint64_t key) {
	// uint64_t ab=abundance_mini[key];
	// uint64_t ab=0;
	// ab<<=32;
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return (uint64_t)key;
}





/** Get the minimizer from a sequence and modify the position parameter.
  * WARNING: If the sequence contains a multiple time the minimizer, return the minimizer position generating a kmer with the longest prefix.
	*
	* @param seq Binariesed sequence
	* @param k size of the kmer
	* @param min_position Position of the minimizer in the kmer (from the end of the suffix). Negative with complement A 1 if the sequence contains multiple values of the minimizer.
	* @param m Minimizer size
	*
	* @return The minimizer value. Negative number if the minimizer is on the reverse complement.
  */
uint64_t get_minimizer(kint seq, const uint8_t k, uint8_t& min_position, const uint8_t m, bool & reversed, bool & multiple) {
	// Init with the first possible minimizer
	uint64_t mini, mmer;
	uint64_t fwd_mini = seq % (1 << (m*2));
	mini = mmer = canonize(fwd_mini, m);
	uint64_t hash_mini = bfc_hash_64(mmer,((1<<(2*m))-1));
	// print_kmer(mini, m); 	cout << endl;
	// cout << hash_mini << endl;

	// Update values regarding the minimizer
	reversed=(mini!=fwd_mini);
	min_position = 0;
	multiple = false;

	// Search in all possible position (from 1) the minimizer
	for (uint8_t i=1; i <= k - m; i++) {
		seq >>= 2;
		fwd_mini = seq % (1 << (m*2));
		mmer = canonize(fwd_mini, m);
		uint64_t hash = bfc_hash_64(mmer,((1<<(2*m))-1));
		// print_kmer(mmer, m); 	cout << endl;
		// cout << (uint)i << " " << hash << endl;

		if (hash_mini > hash) {
			min_position = i;
			mini = mmer;
			reversed=(mini!=fwd_mini);		
			hash_mini = hash;
			multiple = false;
		} else if (hash_mini == hash) {
			min_position = i;
			reversed=(mini!=fwd_mini);

			multiple = true;
		}
	}

	//return hash_mini;
	return mini;
}



char revCompChar(char c) {
	switch (c) {
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
	}
	return 'A';
}



string revComp(const string& s) {
	string rc(s.size(), 0);
	for (int i((int)s.length() - 1); i >= 0; i--) {
		rc[s.size() - 1 - i] = revCompChar(s[i]);
	}
	return rc;
}



string getCanonical(const string& str) {
	return (min(str, revComp(str)));
}



uint64_t nuc2int(char c) {
	return (c >> 1) & 0b11;
}



uint64_t nuc2intrc(char c) {
	return ((c >> 1) & 0b11) ^ 2;
}



inline void updateK(kint& min, uint8_t nuc, kint mask) {
	min <<= 2;
	min += nuc;
	min &= mask;
}



inline void updateM(uint64_t& min, uint8_t nuc, uint64_t mask) {
	min <<= 2;
	min += nuc;
	min &= mask;
}



inline void updateRCK(kint& min, char nuc, uint8_t k) {
	min >>= 2;
	min += ((kint)nuc << (2 * k - 2));
}



inline void updateRCM(uint64_t& min, char nuc, uint8_t m) {
	min >>= 2;
	min += (nuc << (2 * m - 2));
}



/**
  * Read a string for the first k-1 nucleotides and init the kmer and rc_kmer values
  */
void init_kmer(const string & seq, uint64_t & seq_idx, kint & kmer_seq, kint & rc_kmer_seq, const uint8_t k, const kint k_mask) {
	kmer_seq = 0;
	rc_kmer_seq = 0;

	for (uint8_t seq_idx=0 ; seq_idx<k ; seq_idx++) {
		auto nuc = nuc2int(seq[seq_idx]);
		updateK(kmer_seq, nuc, k_mask);
		nuc = nuc ^ 2;
		updateRCK(rc_kmer_seq, nuc, k);
	}
}



void update_kmer(const char nucl, kint & kmer_seq, kint & rc_kmer_seq, const uint8_t k, const 								kint k_mask) {
	auto nuc = nuc2int(nucl);
	updateK(kmer_seq, nuc, k_mask);
	nuc = nuc ^ 2;
	updateRCK(rc_kmer_seq, nuc, k);
}



SuperKmerEnumerator::SuperKmerEnumerator(string & s, const uint8_t k, const uint8_t m)
: seq( s ), seq_idx( 0 )
, k( k ), k_mask( ((kint)1 << (2*k)) - 1 )
, m( m ), m_mask( ((kint)1 << (2*m)) - 1 )
, saved( false ), saved_kmer( kmer_full((kint)0,(uint8_t)0, (uint8_t)0, false) )
, current_kmer( 0 ), current_rc_kmer( 0 )
, mini_candidate( 0 ), rc_mini_candidate( 0 )
, mini_pos ( 1 )
{}



kint SuperKmerEnumerator::next(vector<kmer_full> & kmers) {
	bool to_return = false;
	uint64_t return_val;

	// If start of the sequence, init the kmer and the minimizer
	if (seq_idx == 0) {
		// Init kmer
		init_kmer(seq, seq_idx, current_kmer, current_rc_kmer, k-1, k_mask>>2);
		current_rc_kmer <<= 2;
		auto start = seq_idx - m;
		init_kmer(seq.substr(k-m-1, m), start, mini_candidate, rc_mini_candidate, m, m_mask);

		// Init real minimizer
		mini = get_minimizer(current_kmer, k-1, mini_pos, m, reversed, multiple);
		min_hash = bfc_hash_64(mini,m_mask);
	} 

	if (saved) {
		saved = false;
		saved_kmer.hash_kmer_body(m,m_mask);
		kmers.push_back(move(saved_kmer));
	}

	// Loop over all kmers
	uint64_t line_size = seq.size();
	for (; k-1+seq_idx < line_size; ++seq_idx) {
		// Update the current kmer and the minimizer candidate
		update_kmer(seq[k-1+seq_idx], current_kmer, current_rc_kmer, k, k_mask);
		update_kmer(seq[k-1+seq_idx], mini_candidate, rc_mini_candidate, m, m_mask);
		mini_pos += 1;
		
		// Get canonical minimizer
		kint candidate_canon = min(mini_candidate, rc_mini_candidate);
		uint64_t current_hash = bfc_hash_64((uint64_t)candidate_canon,m_mask);

		//the previous MINIMIZER is outdated
		if (mini_pos > k-m) {
			// Save previous kmers from superkmer
			if (reversed){
				reverse(kmers.begin(), kmers.end());
			}
			to_return = true;
			return_val = min_hash;

			// Prepare new minimizer
			mini = get_minimizer(current_kmer, k, mini_pos, m, reversed, multiple);
			min_hash = bfc_hash_64(mini,m_mask);
		}
		// New minimizer
		else if (current_hash < min_hash) {
			// Save previous kmers from superkmer
			if (reversed){
				reverse(kmers.begin(), kmers.end());
			}
				
			to_return = true;
			return_val = min_hash;

			// Update for the new minimizer value
			mini = candidate_canon;
			mini_pos = 0;
			min_hash = bfc_hash_64(mini,m_mask);
		
			reversed = mini == rc_mini_candidate;
			multiple = false;
		}
		// Equal minimizer
		else if (current_hash == min_hash) {
			// Save previous kmers from superkmer
			if (reversed){
				reverse(kmers.begin(), kmers.end());
			}
				
			to_return = true;
			return_val = min_hash;

			multiple = true;
		}

		// Multiple minimizer is directed regarding the complete kmer value because there is no reference due to the multiple index for minimizers.
		if (multiple) {
			if (current_kmer > current_rc_kmer) {
				reversed = true;
			} else {
				reversed = false;
			}
		}

		if (not reversed) {
			saved_kmer = kmer_full(current_kmer, mini_pos, m, multiple);
		} else {
			saved_kmer = kmer_full(current_rc_kmer, k - m - mini_pos, m, multiple);
		}

		if (to_return and seq_idx > 0) {
			seq_idx += 1;
			saved = true;
			return return_val;
		} else {
			if (seq_idx == 0){
				to_return = false;
			}
			saved_kmer.hash_kmer_body(m,m_mask);
			kmers.push_back(move(saved_kmer));
		}
	}

	if (kmers.size() > 0) {
		if (reversed){
			reverse(kmers.begin(), kmers.end());
		}
		
		return min_hash;
		return mini;
	}
	return (kint)0;
}

