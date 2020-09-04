#include "SuperKmerLight.hpp"
#include "Kmers.hpp"



using namespace std;


/** Construct a superkmer from one kmer and the minimizer position.
	* @param kmer The unsigned int used to represent the binary kmer. The minimizer is not present
	* @param mini_idx The minimizer position in the kmer.
	*/
SKCL::SKCL(kint kmer, const uint8_t mini_idx, const uint32_t indice_v) {
	memset(nucleotides,0,SKCL::allocated_bytes);
	Pow2<kint> anc(2*compacted_size-8);
	for(uint i(0);i<(compacted_size/4);i++){
		nucleotides[SKCL::allocated_bytes-1-i]=kmer/anc;
		kmer%=anc;
		anc>>=8;
	}
	if(compacted_size%4!=0){
		nucleotides[SKCL::allocated_bytes-1-(compacted_size/4)]=(kmer<<(2*(4-compacted_size%4)));
	}
	this->size          = 1;
	indice_value=indice_v;
	this->minimizer_idx = mini_idx;

	this->bytes_used=ceil(static_cast<float>(k - minimizer_size)/4.);
};



uint SKCL::which_byte(uint i){
	//~ return (SKCL::allocated_bytes-i/4);
	return (SKCL::allocated_bytes-1-(ceil((float)i/4)));
}



/**
  * Return the byte index corresponding to the nucletide position.
  * Position 0 is the first nucleotide of the prefix.
  * The minimizer nucleotides doesn't count.
  * 
  * @param position Nucleotide position in the sequence
  *
  * @return The Byte index in the datastructure.
  * The first 4 nucleotides are inside of the last byte of the byte array (little endian style).
  */
uint SKCL::byte_index(uint position){
	return SKCL::allocated_bytes - 1 - position / 4;
}

/**
  * Get the nucleotide value at the position in parameter.
  * Position 0 is the first nucleotide of the prefix.
  * The minimizer nucleotides doesn't count.
  */
// uint8_t SKCL::get_nucleotide(uint8_t position) {
// 	uint byte_pos = byte_index(position);
// 	// cout << byte_pos << endl;
// 	uint8_t nucl = nucleotides[byte_pos];
// 	nucl >>= 2 * (position%4);
// 	nucl &= 0b11;
// 	return nucl;
// }


uint8_t SKCL::get_nucleotide(uint8_t position) {
	cout<<"position:	"<<(int)position<<endl;
	//~ uint8_t compacted_length = k - minimizer_size + size - 1;
	uint8_t byte_pos = allocated_bytes-(position/4)-1;
	cout<<"bp "<<(int)byte_pos<<endl;

	uint8_t nucl = nucleotides[byte_pos];
	print_kmer(nucl,4);
	cout<<endl;
	print_all();
	nucl >>= 2 * (((allocated_bytes - position - 1)%4));
	nucl &= 0b11;
	cout<<endl;
	print_kmer(nucl,1);
	cout<<endl;
	cout<<" end"<<endl;
	return nucl;
}


uint64_t SKCL::interleaved_value(){
	uint64_t value = 0;
	// Suffix interleaved
	uint8_t max_suffix = min((uint)8, (uint)minimizer_idx);
	for (uint8_t i=0 ; i<max_suffix ; i++) {
		uint8_t position = minimizer_idx - 1 - i;
		cout << (uint64_t)i << " " << (uint64_t)position << endl;
		// Get the value of the nucleotide at the position
		uint64_t nucl_value = get_nucleotide(position);
		cout << nucl_value << endl;
		// shift the value to the right place
		nucl_value <<= 62 - i*4;
		// Add the nucleotide to the interleaved
		value |= nucl_value;
	}

	// prefix interleaved
	uint8_t max_prefix = min((uint)8,prefix_size());
	for (uint8_t i=0 ; i<max_prefix ; i++) {
		// Get the nucleotide position
		uint8_t position = minimizer_idx + i;
		// Get the value of the nucleotide at the position
		uint64_t nucl_value = get_nucleotide(position);
		// shift the value to the right place
		nucl_value <<= 60 - i*4;
		// Add the nucleotide to the interleaved
		value |= nucl_value;
	}
	//~ cout<<get_string("	")<<endl;
	//~ print_kmer(value,32);
	//~ cout<<endl;
	//~ cin.get();
	return value;
}



string SKCL::get_string(const string& mini)const {
	string result;
	for(uint i(0);i<bytes_used;++i){
		result+=kmer2str(nucleotides[SKCL::allocated_bytes-1-i],4);
	};
	
	// Remove the nucleotides that are not needed in the last byte
	result=result.substr(0,size+compacted_size-1);
	auto l = result.length();

	// Extract pref-suff (kmer are reversed)
	string suffix(result.substr(l-minimizer_idx, minimizer_idx));
	string prefix(result.substr(0, l-minimizer_idx));
	return(prefix+mini+suffix);
}


kint SKCL::get_ith_kmer(uint ind)const{
	kint result(0);
	
	int skm_nuc(compacted_size+size-1);
	//~ cout<<"ind:	"<<ind<<endl;
	int start (which_byte(compacted_size+ind-1));
	//~ cout<<"start:	"<<start<<endl;
	int length_to_read (which_byte(ind)-start+1);
	//~ cout<<"length_to_read:	"<<length_to_read<<endl;
	//~ int last_bytes(which_byte(ind));
	//~ cout<<"start:	"<<start<endl;
	
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
		cout<<"IDX DEAD"<<endl;
		return false;
	}
	print_kmer( get_ith_kmer(size-start_idx),k);
	cout<<endl;
	print_kmer( kmer.get_compacted(),k);
	cout<<endl;

	return get_ith_kmer(size-start_idx) == kmer.get_compacted();//THE GET COMPACTED SHOULD BE MADE ABOVE
}



uint32_t SKCL::query_kmer_hash(const kmer_full& kmer)const {
	//~ cout<<"query_kmer_hash"<<endl;	
	if(this->query_kmer_bool(kmer)){
		//~ cout<<"query kmer bool true"<<endl;
		return indice_value+kmer.minimizer_idx - (minimizer_idx-size+1);
	}
	//~ cout<<"query kmer bool FAIL"<<endl;
	return -1;
}



bool SKCL::compact_right(const kmer_full& kmf) {
	kint super_kmer_overlap(get_right_overlap());
	kint kmer_overlap(kmf.get_compacted());
	int nuc(kmer_overlap%4);
	kmer_overlap>>=2;
	kmer_overlap%=((kint)1<<(2*(k-1-minimizer_size)));
	if(super_kmer_overlap==kmer_overlap){
		int byte_to_update(SKCL::allocated_bytes-1-((compacted_size+size-1)/4));
		int padding((4-((compacted_size+size)%4))%4);
		nucleotides[byte_to_update] += (nuc<<(2*padding));
		size++;
		// Size of a kmer - size of minimizer + 1 for each supplementary nucleotide
		bytes_used=ceil(static_cast<float>(k - minimizer_size + size - 1)/4.);
		minimizer_idx++;
		return true;
	}	
	return false;
}


uint SKCL::suffix_size()const{
	return (size+compacted_size-minimizer_idx);
}


uint SKCL::prefix_size()const{
	return (minimizer_idx);
}



uint SKCL::interleaved_size()const{
	return min(prefix_size(),suffix_size())*2;
}



bool  SKCL::is_lex_inferior(const SKCL& kmf)const{
	cout<<"is_lex_inferior"<<endl;
	uint64_t interleaved_index(interleaved);
	uint64_t interleaved_query(kmf.interleaved);
	uint interleaved_size_index(interleaved_size());
	uint interleaved_size_query(kmf.interleaved_size());
	cout<<interleaved_size_index<<" "<<interleaved_size_query<<endl;
	//~ if(interleaved_size_index>interleaved_size_query){
		//~ interleaved_index>>=(2*(interleaved_size_index-interleaved_size_query));
	//~ }
	//~ if(interleaved_size_index<interleaved_size_query){
		//~ interleaved_query>>=(2*(interleaved_size_query-interleaved_size_index));
	//~ }
	return true;
	print_kmer(interleaved_index,31);
	cout<<endl;
	print_kmer(interleaved_query,31);
	cout<<endl;
	cin.get();
	return interleaved_index<=interleaved_query;
}




bool  SKCL::suffix_is_prefix(const kmer_full& kmf)const{
	kint suffix_kmer(kmf.suffix);
	kint suffix_superkmer(this->get_suffix());
	print_kmer(suffix_kmer,k);
	cout<<endl;
	print_kmer(suffix_superkmer,k);
	cout<<endl;
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
		cout<<SKCL::allocated_bytes-1-i;
		print_kmer(nucleotides[SKCL::allocated_bytes-1-i],4);cout<<" ";
		
	}
}




