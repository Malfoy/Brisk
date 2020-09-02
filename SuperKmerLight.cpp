#include "SuperKmerLight.hpp"
#include "Kmers.hpp"



using namespace std;



/** Construct a superkmer from one kmer and the minimizer position.
	* @param kmer The unsigned int 64 used to represent the binary kmer.
	* @param mini_idx The minimizer position in the kmer.
	*/
SKCL::SKCL(kint kmer, const uint8_t mini_idx, const uint32_t indice_v) {
	//~ cout<<"CONSTRUCTOR"<<endl;
	memset(nucleotides,0,byte_nuc+1);
	//~ print_kmer(kmer,compacted_size);cout<<endl;
	// print_kmer(kmer,63);cout<<endl;
	Pow2<kint> anc(2*compacted_size-8);
	for(uint i(0);i<(compacted_size/4);i++){
		//~ cout<<"loop:	"<<byte_nuc-i<<endl;
		//~ print_kmer(kmer/anc,4);cout<<endl;
		nucleotides[byte_nuc-i]=kmer/anc;
		//~ print_kmer(nucleotides[byte_nuc-i],4);cout<<endl;
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
	//~ print_all();
	//~ cout<<get_string()<<endl;
};




string SKCL::get_string(const string& mini)const {
	//~ cout<<"GET STRING"<<endl;
	//~ cout<<"lidx du q	"<<(int)minimizer_idx<<endl;
	string result;
	for(uint i(0);i<bytes_used;++i){
		result+=kmer2str(nucleotides[byte_nuc-i],4);
	}
	result+=kmer2str(nucleotides[byte_nuc-bytes_used],4);
	
	//~ cout<<"ALL SEQ GET STRING"<<endl;
	//~ cout<<result<<endl;
	//~ cout<<"wierd add ?:	"<<kmer2str(nucleotides[byte_nuc-bytes_used],4)<<endl;
	//~ cout<<"wierd add ?:	"<<kmer2str(nucleotides[byte_nuc-bytes_used]>>((4-(k+size-1)%4)%4),4)<<endl;
	//~ cout<<"total "<<result<<endl;
	result=result.substr(0,size+compacted_size-1);
	//~ cout<<"total "<<result<<endl;
	string suffix(result.substr(result.size()-minimizer_idx));
	//~ cout<<result.size()-suffix.size()<<endl;
	string prefix(result.substr(0,result.size()-suffix.size()));
	//~ cout<<"suffix:	"<<suffix<<endl;
	//~ cout<<"prefix:	"<<prefix<<endl;
	//~ cin.get();
	return(prefix+mini+suffix);
}


uint which_byte(uint i){
	//~ cout<<"WB "<<i<<" "<<(ceil((float)i/4))<<endl;
	return (byte_nuc+1-(ceil((float)i/4)));
}



kint SKCL::get_ith_kmer(uint ind)const{
	//~ cout<<"GET KMER:	"<<ind<<endl;
	//~ print_all();
	//~ cout<<endl;
	kint result(0);
	
	int skm_nuc(compacted_size+size-1);
	int start (which_byte(compacted_size+ind-1));
	int length_to_read (which_byte(ind)-start+1);
	int last_bytes(which_byte(ind));
	
	//~ cout<<"compacted+ind "<<compacted_size+ind<<endl;
	//~ cout<<"last byte to read"<<last_bytes<<endl;
	//~ cout<<"SKM NUC"<<" "<<skm_nuc<<endl;
	//~ cout<<"length_to_read:	"<<length_to_read<<endl;
	//~ cout<<"start:	"<<start<<endl;
	
	memcpy(&result, &nucleotides[start], length_to_read);
	//~ cout<<(int)result<<endl;
	//~ cout<<"what I read"<<endl;
	//~ print_kmer(result,k);cout<<endl;
	
	int offset( (4-(compacted_size+ind-1)%4)%4);
	//~ cout<<"offset:	"<<offset<<endl;
	result>>=2*offset;
	result&=compact_mask;
	//~ cout<<"get_ith_kmer result"<<endl;
	//~ print_kmer(result,k);cout<<endl;
	//~ cout<<"compacted_size:	"<<compacted_size<<endl;
	//~ print_kmer(result,compacted_size);cout<<endl;
	//~ cout<<"END GET KMER"<<endl;
	//~ cin.get();
	//~ string superkmerseq(get_string(""));
	//~ if(kmer2str(result,compacted_size)!=superkmerseq.substr(ind-1,compacted_size)){
		//~ cout<<"WRONG"<<endl;
		//~ cout<<kmer2str(result,compacted_size)<<"\n"<<superkmerseq.substr(ind-1,compacted_size)<<endl;
		//~ cout<<superkmerseq<<endl;
		//~ cin.get();
	//~ }
	return result;
}



kint SKCL::get_right_overlap()const {
	//~ cout<<"GET RIGHT OVERLAP"<<endl;
	kint result(get_ith_kmer(size));
	//~ print_kmer(result,k-1);cout<<endl;
	result%=((kint)1<<(2*(k-1-minimizer_size)));
	//~ print_kmer(result,k-1);cout<<endl;
	//~ cout<<"END GET RIGHT OVERLAP"<<endl;
	//~ cin.get();
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
	//~ cout<<"query_kmer_bool"<<endl;
	int64_t start_idx  = (int64_t)this->minimizer_idx - (int64_t)kmer.minimizer_idx;
	if(start_idx<0 or (start_idx>=this->size)){
		//~ cout<<"LE FAILT"<<endl;
		return false;
	}
	//~ print_all();
	//~ cout<<"GET kmer number:	"<<(int)size-start_idx<<endl;
	//~ print_kmer(get_ith_kmer(size-start_idx),k);
	//~ cout<<endl;
	//~ print_kmer(kmer.get_compacted(),k);
	//~ cout<<endl;

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
	//~ cout<<"COMPACT RIGHT"<<endl;
	kint super_kmer_overlap(get_right_overlap());
	kint kmer_overlap(kmf.get_compacted());
	//~ print_kmer(kmer_overlap,k);cout<<endl;
	int nuc(kmer_overlap%4);
	//~ cout<<"nuc:	"<<nuc<<endl;
	kmer_overlap>>=2;
	kmer_overlap%=((kint)1<<(2*(k-1-minimizer_size)));
	//~ cout<<"compact_righté"<<endl;

	//~ print_kmer(super_kmer_overlap,k-1-minimizer_size);cout<<endl;
	//~ print_kmer(kmer_overlap,k-1-minimizer_size);cout<<endl;

	if(super_kmer_overlap==kmer_overlap){
		//~ cout<<"YES"<<endl;
		//~ cout<<byte_nuc<<" "<<((compacted_size+size-1))<<endl;
		int byte_to_update(byte_nuc-((compacted_size+size-1)/4));
		//~ cout<<"byte to update:	"<<byte_to_update<<endl;
		int padding((4-((compacted_size+size)%4))%4);
		//~ cout<<"paddin to add:	"<<padding<<endl;
		nucleotides[byte_to_update] += (nuc<<(2*padding));
		size++;
		bytes_used=ceil(((float)size+(float)k)/4);
		minimizer_idx++;
		//~ print_all();
		//~ cout<<"compaction OK, size:	"<<(int)size<<endl;
		//~ cout<<get_string("")<<endl;
		//~ cin.get();
		return true;
	}
	//~ cout<<"NO"<<endl;print_all();
	//~ cout<<"compaction FAILED, size:	"<<(int)size<<endl;
	//~ cin.get();
	
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
	//~ cout<<"PRINT ALL		";
	// for(uint i(0);i<(2*k-minimizer_size)/4;++i){
	// 	print_bin(nucleotides[i]);
	// 	cout<<",";
	// }
	//~ cout<<"seven"<<(int)nucleotides[7]<<endl;
	for(uint i(0);i<(compacted_size+size)/4+1;++i){
		cout<<byte_nuc-i;
		print_kmer(nucleotides[byte_nuc-i],4);cout<<" ";
		//~ cout<<(int)nucleotides[byte_nuc-i]<<endl;
		
	}
	//~ cout<<"\n";
	//~ cout<<endl;
	//~ cout<<"indice_value:	"<<(int)indice_value<<endl;
	//~ cout<<"size:	"<<(int)size<<endl;
	//~ cout<<"minimizer_idx:	"<<(int)minimizer_idx<<endl;
	//~ cout<<"END PRINT ALL"<<endl;
}




