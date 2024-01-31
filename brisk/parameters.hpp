#include "Decycling.h"

#ifndef PARAMS_H
#define PARAMS_H

class Parameters {
public:
	uint8_t k;
	uint8_t m;
	uint8_t b;
	uint8_t m_reduc;
	uint allocated_bytes;
	uint8_t compacted_size;
	uint64_t mask_large_minimizer;
	DecyclingSet* dede;

	/** Datastruct to hold all the parameters needed for Brisk
	 *  @param k kmer size
	 *  @param m minimizer size [1, k].
	 *  @param b The order of magnitude of the number of buckets [1, m]. The number of buckets will be 4^b with 4^(m - bucket_number) minimizers per bucket.
	 **/
	Parameters (uint8_t k, uint8_t m, uint8_t b) {
		this->k = k;
		this->m = m;
		this->b = b;
		this->m_reduc = m-b;
		this->mask_large_minimizer =((uint64_t)1<<(2*m))-1;
		this->compacted_size = k-b;
		this->allocated_bytes = ceil(((double)(2*k-m -b ))/4);
		DecyclingSet* dd = new DecyclingSet(m);
		this->dede = dd;
	}
};

#endif
