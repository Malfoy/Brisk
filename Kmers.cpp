#include <iostream>
#include "Kmers.hpp"
#include "pow2.hpp"



using namespace std;


// RC functions
uint64_t rcbc(uint64_t in, uint64_t n);
// Hash functions (TODO: move them)
uint64_t hash64shift(uint64_t key);



// ----- Kmer class -----
kmer_full::kmer_full(kint value, uint8_t minimizer_idx, uint8_t m, bool multiple_mini) {
	this->minimizer_idx = minimizer_idx;
	this->kmer_s = value;
	this->prefix=(value);
	uint64_t shift((minimizer_idx + m)*2);
	this->prefix>>=shift;
	this->suffix=(value);
	shift=(minimizer_idx)*2;
	this->suffix%=((kint)1<<shift);
	this->multi_mini = multiple_mini;
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



// SUFFIX IS AT RIGHT!!!!!!DO NOT CHANGE THIS
kint kmer_full::get_compacted() const {
	kint result;
	result=prefix;
	result<<=(minimizer_idx*2);
	result+=suffix;
	return result;
}



bool kmer_full::contains_multi_minimizer() const {
	return this->multi_mini;
}



// ----- Useful binary kmer functions -----
template<typename T>
void print_kmer(T num, uint8_t n){
	num &= ((T)1 << (2*n)) - 1;
	T anc((T)1<<(2*(n-1)));
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
uint64_t get_minimizer(kint seq, uint8_t k, uint8_t& min_position, uint8_t m, bool & reversed, bool & multiple) {
	// Init with the first possible minimizer
	uint64_t mini, mmer;
	uint64_t fwd_mini = seq % (1 << (m*2));
	mini = mmer = canonize(fwd_mini, m);
	uint64_t hash_mini = hash64shift(mmer);

	// Update values regarding the minimizer
	reversed=(mini!=fwd_mini);
	min_position = 0;
	multiple = false;

	// Search in all possible position (from 1) the minimizer
	for (uint8_t i=1; i <= k - m; i++) {
		seq >>= 2;
		fwd_mini = seq % (1 << (m*2));
		mmer = canonize(fwd_mini, m);
		uint64_t hash = (hash64shift(mmer));

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

	// cout << endl;
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
void init_kmer(const string & seq, kint & kmer_seq, kint & rc_kmer_seq, const uint8_t k, const kint k_mask) {
	kmer_seq = 0;
	rc_kmer_seq = 0;

	for (uint8_t i=0 ; i<k ; i++) {
		auto nuc = nuc2int(seq[i]);
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


void string_to_kmers_by_minimizer(string & seq, vector<vector<kmer_full> > & kmers, uint8_t k, uint8_t m) {
	// Useful precomputed values
	kint k_mask = ((kint)1 << (2*k + 1)) - 1;
	kint m_mask = ((kint)1 << (2*m + 1)) - 1;

	// Init kmer
	kint current_kmer = 0;
	kint current_rc_kmer = 0;
	init_kmer(seq, current_kmer, current_rc_kmer, k-1, k_mask>>2);
	current_rc_kmer <<= 2;
	// Init minimizer candidates
	kint mini_candidate = 0;
	kint rc_mini_candidate = 0;
	init_kmer(seq.substr(k-m-1, m), mini_candidate, rc_mini_candidate, m, m_mask);

	// Init real minimizer
	bool reversed, multiple;
	uint8_t mini_pos;
	uint64_t mini = get_minimizer(current_kmer, k-1, mini_pos, m, reversed, multiple);
	uint64_t min_hash = hash64shift(mini);

	vector<kmer_full> current_minimizer_kmers = vector<kmer_full>();

	// Loop over all kmers
	uint64_t line_size = seq.size();
	for (uint64_t i = 0; k-1+i < line_size; ++i) {
		// Update the current kmer and the minimizer candidate
		update_kmer(seq[k-1+i], current_kmer, current_rc_kmer, k, k_mask);
		update_kmer(seq[k-1+i], mini_candidate, rc_mini_candidate, m, m_mask);
		// print_kmer(current_kmer, k);
		// cout << endl;

		// Get canonical minimizer
		kint candidate_canon = min(mini_candidate, rc_mini_candidate);
		uint64_t current_hash = hash64shift((uint64_t)candidate_canon);
		

		//the previous MINIMIZER is outdated
		if (mini_pos >= k-m) {
			// cout << "Outdated" << endl;
			// Save previous kmers from superkmer
			if (reversed)
				reverse(current_minimizer_kmers.begin(), current_minimizer_kmers.end());
			kmers.push_back(current_minimizer_kmers);
			current_minimizer_kmers = vector<kmer_full>();

			// Prepare new minimizer
			mini = get_minimizer(current_kmer, k, mini_pos, m, reversed, multiple);
			min_hash = hash64shift(mini);
		}
		// New minimizer
		else if (current_hash < min_hash) {

			// Save previous kmers from superkmer
			if (reversed)
				reverse(current_minimizer_kmers.begin(), current_minimizer_kmers.end());
			kmers.push_back(current_minimizer_kmers);
			current_minimizer_kmers = vector<kmer_full>();

			// Update for the new minimizer value
			mini = candidate_canon;
			mini_pos = 0;
			min_hash = hash64shift(mini);
			reversed = mini == rc_mini_candidate;
			multiple = false;
		}
		// Equal minimizer
		else if (current_hash == min_hash) {
			// Save previous kmers from superkmer
			if (reversed)
				reverse(current_minimizer_kmers.begin(), current_minimizer_kmers.end());
			kmers.push_back(current_minimizer_kmers);
			current_minimizer_kmers = vector<kmer_full>();

			multiple = true;
		}

		mini_pos += 1;
		if (not reversed) {
			current_minimizer_kmers.push_back(kmer_full(current_kmer, mini_pos, m, multiple));
		} else {
			current_minimizer_kmers.push_back(kmer_full(current_rc_kmer, k - m - mini_pos, m, multiple));
		}
		
		// print_kmer(mini, m);
		// cout << endl;
	}

	// Process last kmers
	if (current_minimizer_kmers.size() > 0) {
		if (reversed)
			reverse(current_minimizer_kmers.begin(), current_minimizer_kmers.end());
		kmers.push_back(current_minimizer_kmers);
	}
}
