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
	if(skml.size()-sorted_size>1000000){
		insert_buffer();
	}
	kmers.clear();
}



void Bucket::insert_buffer(){
	for(auto it(skml.begin()+sorted_size);it<skml.end();++it){
		it->interleaved=it->interleaved_value();
	}
	sort(skml.begin()+sorted_size , skml.end(),[ ]( const SKCL& lhs, const SKCL& rhs ){return lhs < rhs;});
	inplace_merge(skml.begin(), skml.begin()+sorted_size ,skml.end() ,[ ]( const SKCL& lhs, const SKCL& rhs ){return lhs < rhs;});
	sorted_size=skml.size();
}



void  Bucket::add_kmers_buffer(vector<kmer_full>& kmers){
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
				//~ cout<<"indice_v from query kmer hash:	"<<(int)indice_v<<endl;
				if (indice_v!=-1) {
					values[indice_v]++;
					++inserted;
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
						break;
					}else{
					}
				}
				if(not isinserted){
					if(skml.size()==skml.capacity()){
						skml.reserve(skml.capacity()*1.5);
					}
					//~ cout<<"LEPUSH	"<<(int)kmer.get_minimizer_idx()<<endl;
					//~ print_kmer(kmer.kmer_s,k);
					//~ cout<<endl;
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
	//OPTIMIZATION POSSIBLE HERE?
	int insert(0);
	for (uint64_t iikk = 0; iikk < kmers.size(); ++iikk) {
		if(find_kmer(kmers[iikk])){
			kmers[iikk].minimizer_idx=69;
			insert++;
		}
	}
	return insert==kmers.size();
}




bool  Bucket::find_kmer_from_interleave(kmer_full& kmer, SKCL& mockskm){
	uint64_t low=lower_bound (skml.begin(), skml.begin()+sorted_size,mockskm,[ ]( const SKCL& lhs, const SKCL& rhs ){return lhs.interleaved < rhs.interleaved;}) - skml.begin();
	while (low<(uint64_t)sorted_size) {
		if(not skml[low].suffix_is_prefix(kmer)){
			break;
		}
		uint32_t indice_v(skml[low].query_kmer_hash(kmer));
		if (indice_v!=-1) {
			values[indice_v]++;
			return true;
		}
	}
	return false;
}



bool Bucket::find_kmer(kmer_full& kmer){
	SKCL mockskm(kmer.kmer_s, kmer.minimizer_idx,0);
	uint size_interleave(min((int)mockskm.size-mockskm.minimizer_idx,(int)mockskm.minimizer_idx)*2);
	if(size_interleave>=4){
		return find_kmer_from_interleave(kmer,mockskm);
	}else{
		if(mockskm.size-mockskm.minimizer_idx>mockskm.minimizer_idx){
			//SUFFIX IS LARGER PREFIX IS MISSING
			if(size_interleave==1){
				for(uint i(0);i<4;++i){
					mockskm.interleaved+=i<<6;
					if(find_kmer_from_interleave(kmer,mockskm)){return true;}
					mockskm.interleaved-=i<<6;
				}
			}
			if(size_interleave==0){
				for(uint ii(0);ii<4;++ii){
					mockskm.interleaved+=ii<<2;
					for(uint i(0);i<4;++i){
						mockskm.interleaved+=i<<6;
						if(find_kmer_from_interleave(kmer,mockskm)){return true;}
						mockskm.interleaved-=i<<6;
					}
					mockskm.interleaved-=ii<<2;
				}
			}
		}else{
			//PREFIX IS LARGER, SUFFIX IS MISSING
			if(size_interleave==1){
				for(uint i(0);i<4;++i){
					mockskm.interleaved+=i<<4;
					if(find_kmer_from_interleave(kmer,mockskm)){return true;}
					mockskm.interleaved-=i<<4;
				}
			}
			if(size_interleave==0){
				for(uint ii(0);ii<4;++ii){
					mockskm.interleaved+=ii;
					for(uint i(0);i<4;++i){
						mockskm.interleaved+=i<<4;
						if(find_kmer_from_interleave(kmer,mockskm)){return true;}
						mockskm.interleaved-=i<<4;
					}
					mockskm.interleaved-=ii;
				}
			}
		}
	}
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
