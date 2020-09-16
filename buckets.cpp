#include <algorithm>
#include "buckets.hpp"
#include "Kmers.hpp"



using namespace std;



void  Bucket::add_kmers(vector<kmer_full>& kmers){
	if(kmers.empty()){return;}
	if(not add_kmers_sorted(kmers)){
		add_kmers_buffer(kmers);
	}
	if(skml.size()-sorted_size>10){
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
				int32_t indice_v(skc.query_kmer_hash(kmer));
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
				if( skml.size()-sorted_size!=0){
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
					int32_t indice_v(skml[i].query_kmer_hash(kmer));
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



bool  Bucket::add_kmers_sorted(vector<kmer_full>& kmers){
	if(sorted_size==0){
		return false;
	}
	uint64_t insert(0);
	for (uint64_t iikk = 0; iikk < kmers.size(); ++iikk) {
		if(kmers[iikk].minimizer_idx!=69){
			insert+=find_kmer(kmers[iikk],kmers);
		}
	}
	return insert==kmers.size();
}




uint  Bucket::find_kmer_from_interleave(kmer_full& kmer, SKCL& mockskm,vector<kmer_full>& kmers){
	uint insertion(0);
	uint64_t low=lower_bound(skml.begin(), skml.begin()+sorted_size,mockskm,[ ]( const SKCL& lhs, const SKCL& rhs ){return lhs.interleaved < rhs.interleaved;}) - skml.begin();
	uint64_t value_max(mockskm.interleaved_value_max());
	while (low<(uint64_t)sorted_size) {
		if(skml[low].interleaved>value_max){
			return 0;
		}
		int32_t indice_v(skml[low].query_kmer_hash(kmer));
		if (indice_v!=-1){
			for(uint i(0);i<kmers.size();++i){
				indice_v=(skml[low].query_kmer_hash(kmers[i]));
				if(indice_v!=-1){
					values[indice_v]++;
					kmers[i].minimizer_idx=69;
					insertion++;
				}
			}
			return insertion;
		}
		low++;
	}
	return 0;
}



uint Bucket::find_kmer(kmer_full& kmer, vector<kmer_full>& v){
	SKCL mockskm(kmer.get_compacted(), kmer.get_minimizer_idx(),0);
	mockskm.interleaved=mockskm.interleaved_value();
	uint prefix_size(mockskm.prefix_size());
	uint suffix_size(mockskm.suffix_size());
	uint size_interleave(min(prefix_size,suffix_size)*2);
	if(size_interleave>=6){
		return find_kmer_from_interleave(kmer,mockskm,v);
	}else{
		if(suffix_size>prefix_size){
			//SUFFIX IS LARGER PREFIX IS MISSING
			if(size_interleave==4){
				for(uint64_t i(0);i<4;++i){
					mockskm.interleaved+=i<<52;
					uint count(find_kmer_from_interleave(kmer,mockskm,v));
					if(count!=0){return count;}
					mockskm.interleaved-=i<<52;
				}
			}
			if(size_interleave==2){
				for(uint64_t ii(0);ii<4;++ii){
					mockskm.interleaved+=ii<<56;
					for(uint64_t i(0);i<4;++i){
						mockskm.interleaved+=i<<52;
						uint count(find_kmer_from_interleave(kmer,mockskm,v));
						if(count!=0){return count;}
						mockskm.interleaved-=i<<52;
					}
					mockskm.interleaved-=ii<<56;
				}
			}
			if(size_interleave==0){
				for(uint64_t iii(0);iii<4;++iii){
					mockskm.interleaved+=iii<<60;
					for(uint64_t ii(0);ii<4;++ii){
						mockskm.interleaved+=ii<<56;
						for(uint64_t i(0);i<4;++i){
							mockskm.interleaved+=i<<52;
							uint count(find_kmer_from_interleave(kmer,mockskm,v));
							if(count!=0){return count;}
							mockskm.interleaved-=i<<52;
						}
						mockskm.interleaved-=ii<<56;
					}
					mockskm.interleaved-=iii<<60;
				}
			}
		}else{
			//PREFIX IS LARGER, SUFFIX IS MISSING
			if(size_interleave==4){
				for(uint64_t i(0);i<4;++i){
					mockskm.interleaved+=i<<54;
					uint count(find_kmer_from_interleave(kmer,mockskm,v));
					if(count!=0){return count;}
					mockskm.interleaved-=i<<54;
				}
			}
			if(size_interleave==2){
				for(uint64_t ii(0);ii<4;++ii){
					mockskm.interleaved+=ii<<58;
					for(uint64_t i(0);i<4;++i){
						mockskm.interleaved+=i<<54;
						uint count(find_kmer_from_interleave(kmer,mockskm,v));
						if(count!=0){return count;}
						mockskm.interleaved-=i<<54;
					}
					mockskm.interleaved-=ii<<58;
				}
			}
			if(size_interleave==0){
				for(uint64_t iii(0);iii<4;++iii){
					mockskm.interleaved+=iii<<62;
					for(uint64_t ii(0);ii<4;++ii){
						mockskm.interleaved+=ii<<58;
						for(uint64_t i(0);i<4;++i){
							mockskm.interleaved+=i<<54;
							uint count(find_kmer_from_interleave(kmer,mockskm,v));
							if(count!=0){return count;}
							mockskm.interleaved-=i<<54;
						}
						mockskm.interleaved-=ii<<58;
					}
					mockskm.interleaved-=iii<<62;
				}
			}
		}
	}
	return 0;
}



void  Bucket::print_kmers(string& result,const  string& mini)const {
	int count(0);
	for(uint64_t isk(0);isk<skml.size();++isk){
		string skm=skml[isk].get_string(mini);
		for (uint64_t i(0); i < skml[isk].size; ++i) {
			result+=skm.substr(i,k)+'	'+to_string(values[skml[isk].indice_value+i])+'\n';
			count+=values[skml[isk].indice_value+i];
			if(check){
				if(real_count[getCanonical(skm.substr(i,k))]!=(int)values[skml[isk].indice_value+i]){
					cout<<"skm:	"<<skm<<endl;
					cout<<(int)values[skml[isk].indice_value+i]<<" "<<i+skml[isk].indice_value<<" "<<(int)skml[isk].size<<endl;
					cout << "minimizer " << mini << endl;
					cout<<skm.substr(i,k)<<" "<<to_string(values[skml[isk].indice_value+i]);
					cout<<"	instead of ";
					cout<<(int)real_count[getCanonical(skm.substr(i,k))]<<endl;
					cout<<"in cursed kmers:	"<<(int)cursed_kmers[str2num(getCanonical((skm.substr(i,k))))/1024][str2num(getCanonical((skm.substr(i,k))))]<<endl;
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



uint64_t Bucket::number_kmer_counted()const{
	uint64_t result(0);
	for(uint64_t i(0);i<values.size();++i){
		result+=values[i];
	}
	return result;
}
