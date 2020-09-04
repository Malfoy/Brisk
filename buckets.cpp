#include <algorithm>
#include "buckets.hpp"
#include "Kmers.hpp"


using namespace std;





void  Bucket::add_kmers(vector<kmer_full>& kmers){
	if(kmers.empty()){return;}
	//~ for(uint i(0);i<kmers.size();++i){
		//~ kmers[i].get_compacted();
	//~ }
	if(not add_kmers_sorted(kmers)){
		add_kmers_buffer(kmers);
	}
	if(skml.size()-sorted_size>000000000){
		insert_buffer();
	}
	kmers.clear();
}



void Bucket::insert_buffer(){
	for(auto it(skml.begin()+sorted_size);it<skml.end();++it){
		it->interleaved=it->interleaved_value();
	}
	sort(skml.begin()+sorted_size , skml.end(),[ ]( const SKCL& lhs, const SKCL& rhs ){return lhs.interleaved < rhs.interleaved;});
	inplace_merge(skml.begin(), skml.begin()+sorted_size ,skml.end() ,[ ]( const SKCL& lhs, const SKCL& rhs ){return lhs.interleaved < rhs.interleaved;});
	sorted_size=skml.size();
}



void  Bucket::add_kmers_buffer(vector<kmer_full>& kmers){
	//~ cout<<"add kmer buffer"<<endl;
	uint64_t inserted(0);
	uint64_t buffsize(skml.size());
	//HERE IF CHECK IF THE KMER ARE IN THE UNSORTED BUFFER
	//FOREACH SUPERKMER
	for (uint64_t i = sorted_size; i <buffsize; i++) {
		SKCL& skc = skml[i];
		//FOREACH KMER
		for (uint64_t ik = 0; ik < kmers.size(); ++ik) {
			kmer_full& kmer = kmers[ik];
			if (kmer.minimizer_idx!=69) {
				uint32_t indice_v(skc.query_kmer_hash(kmer));
				if (indice_v!=-1) {
					values[indice_v]++;
					++inserted;
					cout<<"inserted unsorted"<<endl;
					kmers[ik].minimizer_idx=69;
				}else{
				}
			}
		}
	}
	//HERE WE CREATE NEW SUPERKMERS (OR ellongate THEM)
	if(inserted!=kmers.size()){
		//FOREACH KMER
		for (uint64_t ik = 0; ik < kmers.size(); ++ik) {
			kmer_full& kmer = kmers[ik];
			if (kmer.minimizer_idx!=69) {
				//WE TRY TO COMPACT IT TO THE LAST SUPERKMER
				if(not skml.empty()){
					if(skml[skml.size()-1].compact_right(kmer)){
						if(values.capacity()==values.size()){
							values.reserve(values.size()*1.5);
						}
						cout<<"compact"<<endl;
						values.push_back(1);
						continue;
					}
				}
				bool isinserted=false;
				//FOREACH new SUPERKMER
				for (uint64_t i = buffsize; i < skml.size(); i++) {
					uint32_t indice_v(skml[i].query_kmer_hash(kmer));
					if (indice_v!=-1) {
						values[indice_v]++;
						isinserted=true;
						cout<<"inserted unsorted2"<<endl;
						break;
					}else{
					}
				}
				if(not isinserted){
					if(skml.size()==skml.capacity()){
						skml.reserve(skml.capacity()*1.5);
					}
					cout<<"le push"<<endl;
					skml.push_back(SKCL(kmer.get_compacted(), (int)kmer.get_minimizer_idx(),values.size()));
					if(values.capacity()==values.size()){
						values.reserve(values.size()*1.5);
					}
					values.push_back(1);
				}
			}
		}
	}
}



bool compSKM(const SKCL& s1, const SKCL& s2){
	return s1<s2;
}



//~ bool  Bucket::add_kmers_sorted( vector<kmer_full>& kmers	){
	//~ if(sorted_size==0){
		//~ return false;
	//~ }
	//~ uint64_t inserted(0);
	//~ //FOREACH KMER
	//~ kmer_full kmer = kmers[0];
	//~ SKCL mockskm(kmer.kmer_s, kmer.minimizer_idx,0);
	//~ uint64_t low=lower_bound (skml.begin(), skml.begin()+sorted_size,mockskm,[ ]( const SKCL& lhs, const SKCL& rhs ){return lhs < rhs;}) - skml.begin();
	//~ //FOREACH SUPERKMER
	//~ while (low<(uint64_t)sorted_size) {
		//~ if(not skml[low].suffix_is_prefix(kmer)){
			//~ break;
		//~ }
		//~ //FOREACH KMER
		//~ for (uint64_t iikk = 0; iikk < kmers.size(); ++iikk) {
			//~ if(kmers[iikk].minimizer_idx==69){continue;}
			//~ uint32_t indice_v(skml[low].query_kmer_hash(kmers[iikk]));
			//~ if (indice_v!=-1) {
				//~ values[indice_v]++;
				//~ kmers[iikk].minimizer_idx=69;
				//~ inserted++;
				//~ if(inserted==kmers.size()){
					//~ return true;
				//~ }
			//~ }
		//~ }
		//~ low++;
	//~ }
	//~ return false;
//~ }


bool  Bucket::add_kmers_sorted( vector<kmer_full>& kmers	){
	if(sorted_size==0){
		return false;
	}
	cout<<"add kmer sorted"<<endl;
	//OPTIMIZATION POSSIBLE HERE?
	int insert(0);
	for (uint64_t iikk = 0; iikk < kmers.size(); ++iikk) {
		if(find_kmer(kmers[iikk])){
			kmers[iikk].minimizer_idx=69;
			insert++;
			cout<<"INSERT SORTED"<<endl;
		}else{
			cout<<"FAIL FIND KMER"<<endl;
		}
	}
	return insert==kmers.size();
}




bool  Bucket::find_kmer_from_interleave(kmer_full& kmer, SKCL& mockskm){
	cout<<"find_kmer_from_interleave"<<endl;
	uint64_t low=lower_bound(skml.begin(), skml.begin()+sorted_size,mockskm,[ ]( const SKCL& lhs, const SKCL& rhs ){return lhs.interleaved < rhs.interleaved;}) - skml.begin();
	low=0;
	//~ cout<<low<<endl;
	while (low<(uint64_t)sorted_size) {
		cout<<(int)kmer.minimizer_idx<<endl;
		print_kmer(mockskm.interleaved,32);
		cout<<endl;
		print_kmer(skml[low].interleaved,32);
		cout<<endl;
		
		cin.get();
		//~ if(skml[low].interleaved>mockskm.interleaved){
			//~ cout<<"break suffix is prefix"<<endl;
			//~ break;
		//~ }
		uint32_t indice_v(skml[low].query_kmer_hash(kmer));
		if (indice_v!=-1){
			values[indice_v]++;
			return true;
		}
		low++;
	}
	return false;
}



bool Bucket::find_kmer(kmer_full& kmer){
	cout<<"find_kmer"<<endl;
	cout<<"THIS IS MOCK SM"<<endl;
	
	cout<<"KMER"<<endl;
	print_kmer(kmer.kmer_s,31);
	cout<<endl;
	print_kmer(kmer.get_compacted(),31-8);
	cout<<endl;
	cout<<(int)kmer.minimizer_idx<<endl;
	
	SKCL mockskm(kmer.get_compacted(), kmer.minimizer_idx,0);
	cout<<mockskm.get_string("	mini	")<<endl;
	
	mockskm.interleaved=mockskm.interleaved_value();
	cout<<"interleave"<<endl;
	print_kmer(mockskm.interleaved,32);
	cout<<endl;
	
	uint prefix_size(mockskm.minimizer_idx);
	uint suffix_size((int)mockskm.size+(int)compacted_size-(int)mockskm.minimizer_idx);
	uint size_interleave(min(prefix_size,suffix_size)*2);
	cout<<" size: "<<(int)mockskm.size+(int)compacted_size<<" prefix size"<<prefix_size<<" suffix size"
	<<suffix_size<<endl;
	cout<<"size_interleave "<<size_interleave<<endl;
	if(size_interleave>=0){
		return find_kmer_from_interleave(kmer,mockskm);
	}else{
		if(suffix_size>prefix_size){
			//SUFFIX IS LARGER PREFIX IS MISSING
			if(size_interleave==2){
				for(uint i(0);i<4;++i){
					mockskm.interleaved+=i<<56;
					if(find_kmer_from_interleave(kmer,mockskm)){return true;}
					mockskm.interleaved-=i<<56;
				}
			}
			if(size_interleave==0){
				for(uint ii(0);ii<4;++ii){
					mockskm.interleaved+=ii<<60;
					for(uint i(0);i<4;++i){
						mockskm.interleaved+=i<<56;
						if(find_kmer_from_interleave(kmer,mockskm)){return true;}
						mockskm.interleaved-=i<<56;
					}
					mockskm.interleaved-=ii<<60;
				}
			}
		}else{
			//PREFIX IS LARGER, SUFFIX IS MISSING
			if(size_interleave==2){
				for(uint i(0);i<4;++i){
					mockskm.interleaved+=i<<58;
					if(find_kmer_from_interleave(kmer,mockskm)){return true;}
					mockskm.interleaved-=i<<58;
				}
			}
			if(size_interleave==0){
				for(uint ii(0);ii<4;++ii){
					mockskm.interleaved+=ii<<62;
					for(uint i(0);i<4;++i){
						mockskm.interleaved+=i<<58;
						if(find_kmer_from_interleave(kmer,mockskm)){return true;}
						mockskm.interleaved-=i<<58;
					}
					mockskm.interleaved-=ii<<62;
				}
			}
		}
	}
	cout<<"FAIL"<<endl;
	//~ cin.get();
	return false;
}



void  Bucket::print_kmers(string& result,const  string& mini)const {
	for(uint64_t isk(0);isk<skml.size();++isk){
		string skm=skml[isk].get_string(mini);
		for (uint64_t i(0); i < skml[isk].size; ++i) {
			result+=skm.substr(i,k)+'	'+to_string(values[skml[isk].indice_value+i])+'\n';
			if(check){
				if(real_count[getCanonical(skm.substr(i,k))]!=(int)values[skml[isk].indice_value+i]){
					cout<<"skm:	"<<skm<<endl;
					cout<<(int)values[skml[isk].indice_value+i]<<" "<<i<<" "<<(int)skml[isk].size<<endl;
					cout<<skm.substr(i,k)<<" "<<to_string(values[skml[isk].indice_value+i]);
					cout<<"	instead of ";
					cout<<real_count[getCanonical(skm.substr(i,k))]<<endl;
					counting_errors++;
				}
				real_count[getCanonical(skm.substr(i,k))]=0;
			}
		}
	}
}



uint64_t Bucket::size()const{
	return skml.size();
}



uint64_t Bucket::number_kmer()const{
	uint64_t result(0);
	for(uint64_t i(0);i<skml.size();++i){
		result+=skml[i].size;
	}
	return result;
}
