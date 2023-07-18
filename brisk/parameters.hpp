#ifndef PARAMS_H
#define PARAMS_H

class Parameters {
public:
	uint8_t k;
	uint8_t m;
	uint8_t m_reduc;
	uint8_t m_small;
	uint allocated_bytes;
	uint8_t compacted_size;
	uint64_t mask_large_minimizer;

	/** Datastruct to hold all the parameters needed for Brisk
	 *  @param k kmer size
	 *  @param m minimizer size [1, k].
	 *  @param bucket_magnitude The order of magnitude for bucket number [1, m]. The number of buckets will be 4^bucket_magnitude with 4^(m - bucket_number) minimizer per bucket.
	 **/
	Parameters (uint8_t k, uint8_t m, uint8_t bucket_magnitude) {
		this->k = k;
		this->m = m;
		this->m_reduc = m-bucket_magnitude;
		this->m_small = bucket_magnitude;
		this->mask_large_minimizer =((uint64_t)1<<(2*m))-1;
		this->compacted_size = k - m_small;
		// this->allocated_bytes = 2 * ((compacted_size + 3) / 4);
		this->allocated_bytes = ceil(((double)(2*k-m -m_small ))/4);
	}
};

#endif
