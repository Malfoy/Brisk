#include <iostream>
#include <algorithm>
#include <cmath>
#include "hashing.hpp"
#include "Kmers.hpp"
#include "pow2.hpp"



using namespace std;



// RC functions
uint64_t rcbc(uint64_t in, uint64_t n);
// ----- Kmer class -----
kmer_full::kmer_full(kint value, uint8_t minimizer_idx, uint8_t minimizer_size, DecyclingSet* dede) {
	this->minimizer_idx = minimizer_idx;
	this->kmer_s = value;
	// Shift the kmer to align minizer on the right
	this->minimizer = value >> (2 * minimizer_idx);
	// Mask bits that are not part of the minimizer
	this->minimizer &= ((kint)1 << (2 * minimizer_size))- 1;
	this->dede=dede;
}



kmer_full::kmer_full(){}



kmer_full::kmer_full(kmer_full&& kmer) :
	kmer_s(kmer.kmer_s),
	minimizer(kmer.minimizer),
	minimizer_idx(kmer.minimizer_idx),
	interleaved(move(interleaved)),
	dede(kmer.dede)
{}



kmer_full::kmer_full(const kmer_full&& kmer) :
	kmer_s(kmer.kmer_s),
	minimizer(kmer.minimizer),
	minimizer_idx(kmer.minimizer_idx),
	interleaved(kmer.interleaved),
	dede(kmer.dede)
{}



void kmer_full::copy(const kmer_full& kmer){
	kmer_s=(kmer.kmer_s);
	minimizer=(kmer.minimizer);
	minimizer_idx=(kmer.minimizer_idx);
	dede=(kmer.dede);
	interleaved=(kmer.interleaved);
}



kmer_full & kmer_full::operator=(kmer_full&& kmer) {
	this->minimizer_idx = kmer.minimizer_idx;
	this->kmer_s = kmer.kmer_s;
	this->minimizer = kmer.minimizer;
	this->dede = kmer.dede;
	this->interleaved = move(kmer.interleaved);
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



vector<int> kmer_full::compute_interleaved(const Parameters & params) const {
	vector<int> interleaved;
	uint64_t k = params.k;
	uint64_t pref_reduc = params.m_reduc / 2;
	uint64_t suff_reduc = params.m_reduc - pref_reduc;
	uint64_t suff_size = this->suffix_size() + suff_reduc;
	uint64_t pref_size = k - params.b - suff_size;
	uint64_t max_idx = 2 * max(suff_size, pref_size);
	interleaved.resize(max_idx, -2);
	// Compute suffix interleved
	for (uint offset=0 ; offset<suff_size ; offset++) {
		uint64_t nucl_idx = suff_size - 1 - offset;
		uint8_t byte = ((uint8_t*)&kmer_s)[nucl_idx/4];
		interleaved[2*offset] = lookup[nucl_idx%4][byte];
	}
	// Compute prefix interleved
	for (uint offset=0 ; offset<pref_size ; offset++) {
		uint8_t nucl_idx = suff_size + params.b + offset;
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



kint kmer_full::get_compacted(uint8_t m, uint8_t mini_idx) const {
	kint mask = (((kint)1) << (mini_idx * 2)) - 1;
	kint suffix = kmer_s & mask;
	mask = ~mask;
	kint prefix = (kmer_s >> (2 * m));
	prefix = prefix & mask;
	return prefix + suffix;
}



void kmer_full::replace_slice(kint replacement, size_t position, size_t length){
	// Create a mask for the slice position
	kint mask = (((kint)1) << (2 * length)) - 1;
	replacement &= mask;
	mask <<= 2 * position;
	// Make a hole at the right position
	this->kmer_s &= ~mask;
	// Insert replacement
	replacement <<= 2 * position;
	this->kmer_s += replacement;
}



kint kmer_full::get_unhash_kmer_value(uint8_t m) const {
	// Extract hashed minimizer
	kint minimizer = this->kmer_s >> (this->minimizer_idx * 2);
	kint mask = (((kint)1) << (m * 2)) - 1;
	minimizer &= mask;
	// Unhash mini
	minimizer = bfc_hash_64_inv(minimizer, mask);
	// Recompose value
	mask <<= 2 * this->minimizer_idx;
	minimizer <<= 2 * this->minimizer_idx;
	return (this->kmer_s & ~mask) + minimizer;
}



void kmer_full::unhash_kmer_minimizer(uint8_t m){
	// Extract hashed minimizer
	kint minimizer = this->kmer_s >> (this->minimizer_idx * 2);
	kint mask = (((kint)1) << (m * 2)) - 1;
	minimizer &= mask;
	// Unhash mini
	minimizer = bfc_hash_64_inv(minimizer, mask);
	this->minimizer = minimizer;
	this->replace_slice(minimizer, this->minimizer_idx, m);
}



void kmer_full::hash_kmer_minimizer_inplace(uint8_t m){
	// Extract hashed minimizer
	kint minimizer = this->kmer_s >> (this->minimizer_idx * 2);
	kint mask = (((kint)1) << (m * 2)) - 1;
	minimizer &= mask;
	// Unhash mini
	minimizer = bfc_hash_64(minimizer, mask,dede);
	this->minimizer = minimizer;
	this->replace_slice(minimizer, this->minimizer_idx, m);
}



kmer_full kmer_full::hash_kmer_minimizer_copy(uint8_t m) const {
	kmer_full copy;
	copy.copy(*this);
	copy.hash_kmer_minimizer_inplace(m);
	return copy;
}



// ----- Useful binary kmer functions -----




string kmer2str(__uint128_t num, uint k) {
	if (k == 0)
		return "";
	string res;
	num = num & (((__uint128_t)1 << (2 * k)) - 1);
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



kmer_full * str2kmer(const std::string & str, const uint8_t m,DecyclingSet* dede) {
	kint km_val(str2num(str));
	uint8_t min_pos;
	bool reversed;
	get_minimizer(km_val, str.size(), min_pos, m, reversed, ((kint)1<<(2*m))-1,dede);
	kmer_full * kmer;
	if (not reversed)
		kmer = new kmer_full(km_val, min_pos, m, dede);
	else
		kmer = new kmer_full(km_val, str.size() - m - min_pos, m, dede);
	return kmer;
}



kint str2num2(const string& str) {
	kint res(0);
	for (uint64_t i(1); i <= str.size(); i++) {
		res += ((str[i-1] / 2) % 4)<<(str.size()-i);
	}
	return res;
}



__m128i mm_bitshift_right(__m128i x, unsigned count) {
	//~ assume(count < 128, "count=%u >= 128", count);
	__m128i carry = _mm_srli_si128(x, 8);
	if (count >= 64) return _mm_srli_epi64(carry, count - 64); // the non-carry part is all zero, so return early
	// else
	carry = _mm_slli_epi64(carry, 64 - count);
	x = _mm_srli_epi64(x, count);
	return _mm_or_si128(x, carry);
}


__uint128_t rcb(const __uint128_t& in, uint64_t n) {
	union kmer_u {
		__uint128_t k;
		__m128i m128i;
		uint64_t u64[2];
		uint8_t u8[16];
	};
	kmer_u res = {.k = in};
	static_assert(sizeof(res) == sizeof(__uint128_t), "kmer sizeof mismatch");
	// Swap byte order
	kmer_u shuffidxs = {.u8 = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}};
	_mm_shuffle_epi8(res.m128i, shuffidxs.m128i);
	// Swap nuc order in bytes
	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
	const uint64_t c2 = 0x3333333333333333;
	for (uint64_t& x : res.u64) {
		x = ((x & c1) << 4) | ((x & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
		x = ((x & c2) << 2) | ((x & (c2 << 2)) >> 2); // swap nuc order in 2-nuc
		x ^= 0xaaaaaaaaaaaaaaaa; // Complement;
	}
	// Realign to the right
	res.m128i = mm_bitshift_right(res.m128i, 128 - 2 * n);
	return res.k;
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



kint canonize(kint x, uint64_t n) {
	return min(x, rcb(x, n));
}



bool canonized(kint x, uint64_t n) {
	if (x == canonize(x,n)) {
		return true;
	}
	return false;
}



/** Get the minimizer from a sequence and modify the position parameter.
  * WARNING: If the sequence contains a multiple time the minimizer, return the minimizer position generating a kmer with the longest prefix.
	*
	* @param seq Binariesed sequence
	* @param k size of the kmer
	* @param min_position Position of the minimizer in the kmer (from the end of the suffix). Negative with complement A 1 if the sequence contains multiple values of the minimizer.
	* @param m Minimizer size
	*
	* @return The minimizer value.
  */
uint64_t get_minimizer(kint seq, const uint8_t k, uint8_t& min_position, const uint8_t m, bool & reversed, const uint64_t m_mask,DecyclingSet* dede) {
	// Init with the first possible minimizer
	uint64_t mini, mmer, hash_mini, new_hash;
	uint64_t fwd_mmer = seq & m_mask;
	uint64_t cur_seq = seq;
	mini = canonize(fwd_mmer, m);
	hash_mini = bfc_hash_64(mini,m_mask,dede);
	reversed=(mini!=fwd_mmer);
	min_position = 0;
	// Search in all possible position (from 1) the minimizer
	for (uint64_t i=1; i <= (uint)k - m; i++) {
		cur_seq >>= (kint)2;
		fwd_mmer = ((uint64_t)cur_seq) & m_mask;
		mmer = canonize(fwd_mmer, m);
		new_hash = bfc_hash_64(mmer,m_mask,dede);
		// new mmer is smaller than previous minimizer
		if (new_hash < hash_mini) {
			min_position = i;
			mini = mmer;
			reversed = (mini!=fwd_mmer);
			hash_mini = new_hash;
		// new mmer is equal to previous minimizer
		} else if (hash_mini == new_hash) {
			// new mmer's position is closer to kmer edge than previous minimizer
			if(k-m-i < min_position){
				min_position = k-m-i;
				mini = mmer;
				reversed = (mini!=fwd_mmer);
				hash_mini = new_hash;
			// new mmer's position is as close to kmer edge than previous minimizer
			}else if(k-m-i == min_position){
				// - strand is the canonical strand
				if(!canonized(seq,k)){
					min_position = k-m-i;
					mini = mmer;
					reversed = false;
				}
			}
		}
	}
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



void update_kmer(const char nucl, kint & kmer_seq, kint & rc_kmer_seq, const uint8_t k, const kint k_mask) {
	auto nuc = nuc2int(nucl);
	updateK(kmer_seq, nuc, k_mask);
	nuc = nuc ^ 2;
	updateRCK(rc_kmer_seq, nuc, k);
}



SuperKmerEnumerator::SuperKmerEnumerator(string & s, const uint8_t k, const uint8_t m,DecyclingSet* dd)
: seq( s ), seq_idx( 0 )
, dede(dd)
, k( k ), k_mask( ((kint)1 << (2*k)) - 1 )
, m( m ), m_mask( ((kint)1 << (2*m)) - 1 )
, saved( false ), saved_kmer( kmer_full((kint)0,(uint8_t)0, (uint8_t)0,dd) )
, current_kmer( 0 ), current_rc_kmer( 0 )
, mini_candidate( 0 ), rc_mini_candidate( 0 )
, mini_pos ( 1 ), mini_hash ( 0 )
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
		mini = get_minimizer(current_kmer, k-1, mini_pos, m, reversed, m_mask,dede);
		mini_hash = bfc_hash_64((uint64_t)mini,m_mask,dede);
	} 
	if (saved) {
		saved = false;
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
		uint64_t current_hash = bfc_hash_64((uint64_t)candidate_canon,m_mask,dede);
		//the previous MINIMIZER is outdated
		if (mini_pos > k-m) {
			// Save previous kmers from superkmer
			if (reversed){
				reverse(kmers.begin(), kmers.end());
			}
			to_return = true;
			return_val = mini;
			// Prepare new minimizer
			mini = get_minimizer(current_kmer, k, mini_pos, m, reversed, m_mask,dede);
			mini_hash = bfc_hash_64((uint64_t)mini,m_mask,dede);
		}
		// New minimizer
		else if (current_hash < mini_hash) {
			// Save previous kmers from superkmer
			if (reversed){
				reverse(kmers.begin(), kmers.end());
			}
			to_return = true;
			return_val = mini;
			// Update for the new minimizer value
			mini_hash = current_hash;
			mini_pos = 0;
			reversed = (candidate_canon == rc_mini_candidate);
		}
		if (not reversed) {
			saved_kmer = kmer_full(current_kmer, mini_pos, m, dede);
		} else {
			saved_kmer = kmer_full(current_rc_kmer, k - m - mini_pos, m, dede);
			saved_kmer.minimizer = mini;
		}
		if (to_return and seq_idx > 0) {
			seq_idx += 1;
			saved = true;
			return return_val;
		} else {
			if (seq_idx == 0){
				to_return = false;
			}
			kmers.push_back(move(saved_kmer));
		}
	}
	if (kmers.size() > 0) {
		if (reversed){
			reverse(kmers.begin(), kmers.end());
		}		
		return mini;
	}
	return (kint)0;
}



void SuperKmerEnumerator::update(uint8_t new_m, DecyclingSet* new_dede){
	m = new_m;
	m_mask = ((kint)1 << (2*m)) - 1;
	dede = new_dede;
	saved = false;
	mini_pos = k-m+1;
}
