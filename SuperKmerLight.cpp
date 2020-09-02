#include "SuperKmerLight.hpp"
#include "Kmers.hpp"



using namespace std;



/** Construct a superkmer from one kmer and the minimizer position.
	* @param kmer The unsigned int 64 used to represent the binary kmer.
	* @param mini_idx The minimizer position in the kmer.
	*/
SKCL::SKCL(kint kmer, const uint8_t mini_idx, const uint32_t indice_v) {
	memset(nucleotides,0,byte_nuc+1);
	Pow2<kint> anc(2*compacted_size-8);
	for(uint i(0);i<(compacted_size/4);i++){
		nucleotides[byte_nuc-i]=kmer/anc;
		kmer%=anc;
		anc>>=8;
	}
	if(compacted_size%4!=0){
		nucleotides[byte_nuc-(compacted_size/4)]=(kmer<<(2*(4-compacted_size%4)));
	}
	this->size          = 1;
	indice_value=indice_v;
	this->minimizer_idx = mini_idx;
	this->bytes_used=ceil((float)k/4);
};




string SKCL::get_string(const string& mini)const {
	string result;
	for(uint i(0);i<bytes_used;++i){
		result+=kmer2str(nucleotides[byte_nuc-i],4);
	}
	result+=kmer2str(nucleotides[byte_nuc-bytes_used],4);
	
	result=result.substr(0,size+compacted_size-1);
	string suffix(result.substr(result.size()-minimizer_idx));
	string prefix(result.substr(0,result.size()-suffix.size()));
	return(prefix+mini+suffix);
}


uint which_byte(uint i){
	return (byte_nuc+1-(ceil((float)i/4)));
}



kint SKCL::get_ith_kmer(uint ind)const{
	kint result(0);
	
	int skm_nuc(compacted_size+size-1);
	int start (which_byte(compacted_size+ind-1));
	int length_to_read (which_byte(ind)-start+1);
	int last_bytes(which_byte(ind));
	
	
	memcpy(&result, &nucleotides[start], length_to_read);
	
	int offset( (4-(compacted_size+ind-1)%4)%4);
	result>>=2*offset;
	result&=compact_mask;
	return result;
}



kint SKCL::get_right_overlap()const {
	kint result(get_ith_kmer(size));
	result%=((kint)1<<(2*(k-1-minimizer_size)));
	return result;
}



kint SKCL::get_suffix()const {
	kint result(get_ith_kmer(size));
	result%=((kint)1<<(2*(minimizer_idx)));
	return result;
}




bool  SKCL::operator < (const  SKCL& str) const {
	kint s1((this->get_suffix()));
	kint s2((str.get_suffix()));
	if(minimizer_idx>=str.minimizer_idx){
		s1>>=(2*(minimizer_idx-str.minimizer_idx));
	}else{
		s2>>=(2*(str.minimizer_idx-minimizer_idx));
	}
	if(s1==s2){
		return minimizer_idx < str.minimizer_idx;
	}
	return  s1<s2;
}



bool SKCL::query_kmer_bool(const kmer_full& kmer)const {
	int64_t start_idx  = (int64_t)this->minimizer_idx - (int64_t)kmer.minimizer_idx;
	if(start_idx<0 or (start_idx>=this->size)){
		return false;
	}
	return get_ith_kmer(size-start_idx) == kmer.get_compacted();//THE GET COMPACTED SHOULE BE MADE ABOVE
}



uint32_t SKCL::query_kmer_hash(const kmer_full& kmer)const {
	//~ cout<<"query_kmer_hash"<<endl;
	if(this->query_kmer_bool(kmer)){
		return indice_value+kmer.minimizer_idx - (minimizer_idx-size+1);
	}
	return -1;
}



bool SKCL::compact_right(const kmer_full& kmf) {
	kint super_kmer_overlap(get_right_overlap());
	kint kmer_overlap(kmf.get_compacted());
	int nuc(kmer_overlap%4);
	kmer_overlap>>=2;
	kmer_overlap%=((kint)1<<(2*(k-1-minimizer_size)));
	if(super_kmer_overlap==kmer_overlap){
		int byte_to_update(byte_nuc-((compacted_size+size-1)/4));
		int padding((4-((compacted_size+size)%4))%4);
		nucleotides[byte_to_update] += (nuc<<(2*padding));
		size++;
		bytes_used=ceil(((float)size+(float)k)/4);
		minimizer_idx++;
		return true;
	}	
	return false;
}



bool  SKCL::suffix_is_prefix(const kmer_full& kmf)const{
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



bool  SKCL::suffix_is_prefix(const SKCL& kmf)const{
	kint suffix_kmer(kmf.get_suffix());
	kint suffix_superkmer(this->get_suffix());
	if(minimizer_idx>=kmf.minimizer_idx){
		suffix_superkmer>>=(2*(minimizer_idx-kmf.minimizer_idx));
	}else{
		return false;
		suffix_kmer>>=(2*(kmf.minimizer_idx-minimizer_idx));
	}
	return (suffix_superkmer==suffix_kmer);
}



//DEBUG FUNCTIONS
void print_bin(uint8_t n) {
	uint8_t mask = 1;
	mask <<= 7;
	for (uint64_t i(0); i < 8; ++i) {
		cout << n / mask;
		if (n / mask == 1) {
			n -= mask;
		}
		mask >>= 1;
	}
	// cout << "\n";
}



void SKCL::print_all()const{
	for(uint i(0);i<(compacted_size+size)/4+1;++i){
		cout<<byte_nuc-i;
		print_kmer(nucleotides[byte_nuc-i],4);cout<<" ";
		
	}
}


uint64_t SKCL::interleaved_value()const{
	return 0;
}




