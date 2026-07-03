#pragma once

class flat_chain_t;
#include "flat_params.h"

class hitdata
	{
public:
// Set by reseeker
	const flat_chain_t *query = 0;
	const flat_chain_t *target = 0;
	const char *path = 0;
	uint ncol = 0;
	float nu_self_score = 0;
	float mega_self_score = 0;
	float nu_fwd_score = 0;
	float nu_rev_score = 0;
	float mega_fwd_score = 0;
	float mega_rev_score = 0;
	uint qlo = 0;
	uint tlo = 0;
	float lddt = 0;
	float dali = 0;
	float TS = 0;

// Derived, set by fill()
	uint ids = UINT_MAX;
	uint diffs = UINT_MAX;
	uint gaps = UINT_MAX;
	uint qhi = UINT_MAX;
	uint thi = UINT_MAX;
	double pvalue = FLT_MAX;

private:
	static const uint CIGAR_BUFSIZE = 4000;
	char cigar_buf[CIGAR_BUFSIZE];
	char *cigar_heap = 0;
	uint cigar_heap_cap = 0;
	uint cigar_len = 0;

public:
	~hitdata()
		{
		cigar_free();
		}

	const char *cigar_ptr() const
		{
		return cigar_heap ? cigar_heap : cigar_buf;
		}

	uint cigar_length() const
		{
		return cigar_len;
		}

	void reset()
		{
		query = 0;
		target = 0;
		path = 0;
		ncol = 0;
		nu_self_score = 0;
		mega_self_score = 0;
		nu_fwd_score = 0;
		nu_rev_score = 0;
		mega_fwd_score = 0;
		mega_rev_score = 0;
		lddt = 0;
		dali = 0;
		TS = 0;

		qhi = UINT_MAX;
		thi = UINT_MAX;
		pvalue = FLT_MAX;
		cigar_free();
		cigar_len = 0;
		}

	void fill(const flat_params &params);
	double calc_pvalue(double TS, PVALUE_MODE pvm);

	void cigar_free();
	bool cigar_ensure(uint need);
	void cigar_put_uint(uint n);
	void cigar_put_op(uint n, char op);
	};
