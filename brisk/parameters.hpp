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

	Parameters (uint8_t k, uint8_t m, uint8_t bucket_magnitude) {
		this->k = k;
		this->m = m;
		this->m_reduc = bucket_magnitude;
		this->m_small = m - m_reduc;

		this->compacted_size = k - m_small;
		this->allocated_bytes = ceil((float)(2 * compacted_size) / 4.);
	}
};

#endif