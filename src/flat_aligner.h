#pragma once

class flat_features;

class flat_aligner
	{
public:
	flat_features *m_ff = 0;
	float *__restrict m_pssmT = 0;
	float m_rev_weight = 0;

	string m_labelQ;
	string m_labelT;

	const uint8_t *__restrict m_profQ = 0;
	const uint8_t *__restrict m_profT = 0;

	uint m_LQ = 0;
	uint m_LT = 0;

	const float **__restrict m_scratch_pssms = 0;
	float *m_scratch_rows = 0;
	uint8_t *__restrict m_TB = 0;

	float m_open = -3;
	float m_ext = -1;

	uint m_maxL = 4000;

	float m_score = 0;
	uint m_loQ = UINT_MAX;
	uint m_loT = UINT_MAX;
	char *m_path_buffer = 0;
	uint m_ncol = 0;

public:
	void alloc();
	void freemem();
	void cacheT(const string &labelT, const uint8_t *profT, uint LT);
	void cacheT_reversed(const string &labelT, const uint8_t *profT, uint LT);
	void alignQ(const string &labelQ, const uint8_t *profQ, uint LQ);
	void write_aln(FILE *f) const;
	void write_tsv(FILE *f) const;
	uint get_path_str(string &path) const;
	float get_self_rev_score(const string &labelQ, const uint8_t *profQ, uint LQ);
	};
