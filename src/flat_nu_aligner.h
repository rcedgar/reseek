#pragma once

#pragma once

#include "parasail.h"

class flat_alphas;

class flat_nu_aligner
	{
public:
	static parasail_matrix_t m_matrix;
	static int m_open;
	static int m_ext;
	static int m_saturated_score;

	static atomic<uint32_t> m_aln_count;
	static atomic<uint32_t> m_saturated_count;

public:
	string m_labelQ;
	const uint8_t *m_codeseqQ = 0;
	uint m_LQ = 0;

	string m_labelT;
	const uint8_t *m_codeseqT = 0;
	uint m_LT = 0;

	parasail_profile_t *m_parasail_profQ = 0;
	parasail_result_t *m_result = 0;
	int m_score = 0;
	uint m_loQ = UINT_MAX;
	uint m_loT = UINT_MAX;
	string m_path; // WARNING -- semi-global(?)

public:
	void clear_aln()
		{
		m_score = 0;
		}

	void set_query(
		const string &labelQ,
		const byte *codeseqQ,
		uint LQ);

	int align_score_only(
		const string &labelT,
		const byte *codeseqT,
		uint LT);

	int align_path(
		const string &labelT,
		const byte *codeseqT,
		uint LT);

public:
	static bool init();
	};
