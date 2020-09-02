#include "SuperKmerCount.hpp"
#include "Kmers.hpp"



uint64_t nb_kmer_read(0);



using namespace std;



/** Construct a superkmer from one kmer and the minimizer position.
	* @param kmer The unsigned int 64 used to represent the binary kmer.
	* @param mini_idx The minimizer position in the kmer.
	*/
SKC::SKC(const kint kmer, const uint8_t mini_idx, const uint32_t indice_v) {
	this->sk            = kmer;
	this->size          = 1;
	indice_value=indice_v;
	this->minimizer_idx = mini_idx;
};



//THE PREFIX IS AT RIGHT!
 skint SKC::get_prefix()  const {
	return sk>>(2*minimizer_idx);
}



//THE SUFFIX IS AT LEFT!
 skint SKC::get_suffix() const  {
	return sk%((kint)1<<(2*minimizer_idx));
}



bool SKC::compact_right(const kmer_full& kmf) {
	if ((k+(int)this->size-1-(int)minimizer_size )>= sizeof(skint)*4) {
		return false;
	}
	if(k-1<minimizer_size+minimizer_idx){
		return false;
	}
	skint prefix(get_prefix());
	skint suffix(get_suffix());
	if( (skint) (kmf.suffix>>2)==suffix){
		// if(true	) {
		// cout<<(2*(k-minimizer_size-minimizer_idx-1))<<endl;
		if((skint)kmf.prefix==(prefix%((skint)1<<(2*(k-minimizer_size-minimizer_idx-1))))) {
			sk<<=2;
			sk += (kmf.suffix % 4);
			size++;
			minimizer_idx++;
			return true;
		}else{
		}
	}else{
	}
	return false;
}



/** Used to compact a new nucleotide from a kmer on the right of the superkmer.
	* @param kmer The binary representation of the kmer to compact on the right. Must be on the same strand than the superkmer.
	* @return True if the kmer is inserted false otherwise.
	*/
	//NOT USED
// bool SKC::compact_left(const kint kmer_val) {
// 	kint begin_sk = this->sk >> (this->size * 2);
// 	begin_sk &= (((kint)1 << (2 * (k - 1))) - 1);
// 	kint end_kmer = kmer_val & (((kint)1 << (2 * (k - 1))) - 1);
// 	if (begin_sk == end_kmer) {
// 		this->sk += ((kint)(kmer_val >> (2 * (k - 1)))) << (2 * (k + this->size - 1));
// 		//           Select 2 left bits         Shift to the beginning of the sk
// 		this->size += 1;
// 		return true;
// 	}
// 	return false;
// }



/** Look for the presence of the kmer inside of the superkmer.
	* The function supposes that the minimizer in the kmer is in the same strand than in sk.
	* @param kmer_val The binary value for the kmer.
	* @param minimizer_idx The index of the minimizer in the kmer_val
	* @return True if the kmer is present inside of the sk.
	*/
bool SKC::is_present(kint kmer_val, uint64_t kmer_minimizer_idx) {
	int64_t start_idx  = (int64_t)this->minimizer_idx - (int64_t)kmer_minimizer_idx;
	if(start_idx<0 or (start_idx>=this->size)){
		return false;
	}
	kint aligned_sk = (kint)((this->sk >> (2 * start_idx)) & k_mask);
	return aligned_sk == kmer_val;
}



bool SKC::is_present(kmer_full kmf) {
	if(minimizer_idx>=kmf.minimizer_idx){
		if((int)kmf.minimizer_idx>=(int)minimizer_idx+1-size){
			if(kmf.prefix==get_prefix()%((kint)1<<(2*(k-minimizer_size-kmf.minimizer_idx)))) {
				if(kmf.suffix==get_suffix()>>(2*(minimizer_idx -kmf.minimizer_idx))) {
					return true;
				}
			}
		}
	}
	return false;
}



bool SKC::is_present_brutforce(kmer_full kmer, uint8_t & mini_k_idx) {
	for (mini_k_idx=0 ; mini_k_idx<=k-minimizer_size ; mini_k_idx++) {
		if (this->is_present(kmer.kmer_s, mini_k_idx))
			return true;
	}
	return false;
}



bool SKC::query_kmer_bool(const kmer_full& kmer) {
	if(this->is_present(kmer)){
		return true;
	}
	return false;
}


uint32_t SKC::query_kmer_hash(const kmer_full& kmer) {
	if(this->is_present(kmer)){
		return indice_value+kmer.minimizer_idx - (minimizer_idx-size+1);
	}
	return -1;
}



// --- Pretty printing functions ---
void _out_kmer(ostream& out, kint kmer, uint64_t size) {
	kmer &= (((kint)1) << (2 * size)) - 1;
	kint anc((kint)1 << (2 * (size - 1)));
	for (uint64_t i(0); i < size and anc != 0; ++i) {
		uint64_t nuc = kmer / anc;
		kmer         = kmer % anc;
		if (nuc == 2) {
			out << "T";
		}
		if (nuc == 3) {
			out << "G";
		}
		if (nuc == 1) {
			out << "C";
		}
		if (nuc == 0) {
			out << "A";
		}
		if (nuc >= 4) {
			out << nuc << endl;
			out << "WTF" << endl;
		}
		anc >>= 2;
	}
}



// ostream& operator<<(ostream& out, const SKC& skc) {
// 	// Print the superkmer
// 	_out_kmer(out, skc.sk, k + skc.size - 1);
// 	out << endl;
// 	// Print the minimizer
// 	uint64_t left_spaces = k + skc.size - minimizer_size - skc.minimizer_idx - 1;
// 	for (uint64_t i = 0; i < left_spaces; i++)
// 		out << ' ';
// 	_out_kmer(out, skc.sk >> (2 * skc.minimizer_idx), minimizer_size);
// 	out << " (mini idx " << (uint64_t)skc.minimizer_idx << ")" << endl;
// 	return out;
// }







bool  SKC::suffix_is_prefix(const kmer_full& kmf){
	kint suffix_kmer(kmf.suffix);
	kint suffix_superkmer(this->get_suffix());
	if(minimizer_idx>=kmf.minimizer_idx){
		suffix_superkmer>>=(2*(minimizer_idx-kmf.minimizer_idx));
	}else{
		return false;
		suffix_kmer>>=(2*(kmf.minimizer_idx-minimizer_idx));
	}
	return (suffix_superkmer==suffix_kmer);
}
