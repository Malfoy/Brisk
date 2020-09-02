#include <iostream>
#include "Kmers.hpp"
#include "pow2.hpp"



using namespace std;



const uint64_t k = 63;
const uint64_t minimizer_size = 12;
const uint64_t compacted_size = k-minimizer_size;
const uint64_t super_minimizer_size(minimizer_size+4);
const uint64_t byte_nuc=ceil((2*k-minimizer_size)/4)-1;
uint64_t counting_errors=0;
bool check=false;
robin_hood::unordered_flat_map<string, uint64_t> real_count;
const kint k_mask = (((kint)1) << (2*k)) - 1;
const kint compact_mask = (((kint)1) << (2*compacted_size)) - 1;



string intToString(uint64_t n) {
	if (n < 1000) {
		return to_string(n);
	}
	string end(to_string(n % 1000));
	if (end.size() == 3) {
		return intToString(n / 1000) + "," + end;
	}
	if (end.size() == 2) {
		return intToString(n / 1000) + ",0" + end;
	}
	return intToString(n / 1000) + ",00" + end;
}



// ----- Kmer class -----
uint8_t kmer_full::get_minimizer_idx() const {
	if (this->minimizer_idx < 0) {
		return (uint8_t)(-this->minimizer_idx - 1);
	} else {
		return (uint8_t)this->minimizer_idx;
	}
}



kmer_full::kmer_full(int8_t minimizer_idx, kint value) {
	this->minimizer_idx = minimizer_idx;
	this->kmer_s = value;
	this->prefix=(value);
	uint64_t shift((minimizer_idx+minimizer_size)*2);
	this->prefix>>=shift;
	this->suffix=(value);
	shift=(minimizer_idx)*2;
	this->suffix%=((kint)1<<shift);
}



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



//~ string kmer2str(uint256_t num, uint k) {
	//~ string res;
	//~ Pow2<uint256_t> anc(2 * (k - 1));
	//~ for (uint64_t i(0); i < k; ++i) {
		//~ auto nuc = num / anc;
		//~ num             = num % anc;
		//~ if (nuc == 3) {
			//~ res += "G";
		//~ }
		//~ if (nuc == 2) {
			//~ res += "T";
		//~ }
		//~ if (nuc == 1) {
			//~ res += "C";
		//~ }
		//~ if (nuc == 0) {
			//~ res += "A";
		//~ }
		//~ if (nuc >= 4) {
			//~ cout << "WTF kmer2str" << endl;
			//~ // cout<<(uint6)anc.value()<<endl;
			//~ cout<<nuc<<endl;
			//~ return "";
		//~ }
		//~ anc >>= 2;
	//~ }
	//~ return res;
//~ }



// SUFFIX IS AT RIGHT!!!!!!DO NOT CHANGE THIS
kint kmer_full::get_compacted() const {
	kint result;
	result=prefix;
	result<<=(minimizer_idx*2);
	result+=suffix;
	return result;
}



bool kmer_full::contains_multi_minimizer() const {
	return this->minimizer_idx < 0;
}



// ----- Useful binary kmer functions -----
void print_kmer(__uint128_t num,uint64_t n){
	__uint128_t anc((__uint128_t)1<<(2*(n-1)));
	for(uint64_t i(0);i<n and anc!=0;++i){
		uint64_t nuc=num/anc;
		num=num%anc;
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
	// cout<<endl;
}



kint str2num(const string& str) {
	kint res(0);
	for (uint64_t i(0); i < str.size(); i++) {
		res <<= 2;
		res += (str[i] / 2) % 4;
	}
	return res;
}



kint str2num2(const string& str) {
	kint res(0);
	for (uint64_t i(1); i <= str.size(); i++) {
		res += ((str[i-1] / 2) % 4)<<(str.size()-i);
	}
	return res;
}



__uint128_t rcb(const __uint128_t& in) {
	assume(k <= 64, "k=%u > 64", k);
	union kmer_u {
		__uint128_t k;
		__m128i m128i;
		uint64_t u64[2];
		uint8_t u8[16];
	};
	kmer_u res = {.k = in};
	// static_assert(sizeof(res) == sizeof(__uint128_t), "kmer sizeof mismatch");
	// Swap byte order
	kmer_u shuffidxs = {.u8 = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}};
	res.m128i = _mm_shuffle_epi8(res.m128i, shuffidxs.m128i);
	// Swap nuc order in bytes
	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
	const uint64_t c2 = 0x3333333333333333;
	for (uint64_t& x : res.u64) {
		x = ((x & c1) << 4) | ((x & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
		x = ((x & c2) << 2) | ((x & (c2 << 2)) >> 2); // swap nuc order in 2-nuc
		x ^= 0xaaaaaaaaaaaaaaaa; // Complement;
	}
	// Realign to the right
	res.m128i = mm_bitshift_right(res.m128i, 128 - 2 * k);
	return res.k;
}



uint64_t rcb(const uint64_t& in) {
	assume(k <= 32, "k=%u > 32", k);
	// Complement, swap byte order
	uint64_t res = __builtin_bswap64(in ^ 0xaaaaaaaaaaaaaaaa);
	// Swap nuc order in bytes
	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
	const uint64_t c2 = 0x3333333333333333;
	res               = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
	res               = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc
	// Realign to the right
	res >>= 64 - 2 * k;
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



// It's quite complex to bitshift mmx register without an immediate (constant) count
// See: https://stackoverflow.com/questions/34478328/the-best-way-to-shift-a-m128i
__m128i mm_bitshift_left(__m128i x, unsigned count) {
	assume(count < 128, "count=%u >= 128", count);
	__m128i carry = _mm_slli_si128(x, 8);
	if (count >= 64)                           //TODO: bench: Might be faster to skip this fast-path branch
		return _mm_slli_epi64(carry, count - 64); // the non-carry part is all zero, so return early
	// else
	carry = _mm_srli_epi64(carry, 64 - count);
	x = _mm_slli_epi64(x, count);
	return _mm_or_si128(x, carry);
}



__m128i mm_bitshift_right(__m128i x, unsigned count) {
	assume(count < 128, "count=%u >= 128", count);
	__m128i carry = _mm_srli_si128(x, 8);
	if (count >= 64)
		return _mm_srli_epi64(carry, count - 64); // the non-carry part is all zero, so return early
	// else
	carry = _mm_slli_epi64(carry, 64 - count);
	x = _mm_srli_epi64(x, count);
	return _mm_or_si128(x, carry);
}



uint64_t reversebits(uint64_t b){
	return (((b * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32);
}



Pow2<uint64_t> minimizer_number(2 * super_minimizer_size);



/** Get the minimizer from a sequence and modify the position parameter.*/
int64_t get_minimizer(kint seq, int8_t& min_position) {
	// Init with the first possible minimizer
	int64_t mini, mmer;
	int64_t fwd_mini = seq % minimizer_number;
	mini = mmer = canonize(fwd_mini, super_minimizer_size);
	bool multiple_mini = false;
	bool reversed=(mini!=fwd_mini);
	uint64_t hash_mini = hash64shift(mmer);
	min_position = 0;
	// Search in all possible position (from 1) the minimizer
	for (uint64_t i=1; i <= k - super_minimizer_size; i++) {
		seq >>= 2;
		fwd_mini = seq % minimizer_number;
		mmer = canonize(fwd_mini, super_minimizer_size);
		uint64_t hash = (hash64shift(mmer));
		// cout<<hash<<" "<<hash_mini<<endl;
		if (hash_mini > hash) {
			// cout<<"new pos mini	"<<i<<endl;
			min_position = i;
			mini = mmer;
			reversed=(mini!=fwd_mini);
			hash_mini = hash;
			multiple_mini = false;
		} else if ((hash_mini == hash) and (not multiple_mini)) {
			multiple_mini = true;
			min_position = - min_position - 1;
		}
	}
	if(reversed){
		mini*=-1;
	}
	return ((int64_t)mini);
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
