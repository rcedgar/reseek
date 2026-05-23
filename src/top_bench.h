#pragma once

#include "lookup.h"

class top_bench
	{
public:
	bool m_scores_are_evalues = false;
	float m_topsum3 = FLT_MAX;
	float m_top_SEPQ0_001 = FLT_MAX;
	float m_top_SEPQ0_01 = FLT_MAX;
	float m_top_SEPQ0_1 = FLT_MAX;
	lookup *m_look = 0;
	float *m_score_top_tp = 0;
	float *m_score_top_fp = 0;
	uint *m_domidx_top_tp = 0;
	uint *m_domidx_top_fp = 0;
	float *m_scores = 0;
	bool *m_tps = 0;
	uint *m_order = 0;

public:
	top_bench()
		{
		}

	~top_bench()
		{
		myfree(m_score_top_tp);
		myfree(m_score_top_fp);
		}

public:
	bool better(float score1, float score2) const
		{
		if (m_scores_are_evalues)
			return score1 < score2;
		else
			return score1 > score2;
		}

	void alloc();
	void clear_hits_and_results();
	double bench(const string &msg = "");

	void read_lookup(const string &fn);

	void read_hits(
		const string &fn,
		uint qidx,
		uint tidx,
		uint scoreidx,
		bool triangle);

	void read_tophits(const string &fn);

	void write_top_hits(const string &fn) const;

	bool is_ignored(uint domIdx_i, uint domIdx_j) const
		{
		assert(m_look);
		return m_look->is_ignored_ij(domIdx_i, domIdx_j);
		}

	bool is_tp(uint domIdx_i, uint domIdx_j) const
		{
		assert(m_look);
		if (m_look->is_ignored_ij(domIdx_i, domIdx_j))
			return false;
		return m_look->is_tp_ij(domIdx_i, domIdx_j);
		}

	bool is_fp(uint domIdx_i, uint domIdx_j) const
		{
		assert(m_look);
		if (m_look->is_ignored_ij(domIdx_i, domIdx_j))
			return false;
		return !m_look->is_tp_ij(domIdx_i, domIdx_j);
		}

	float get_missing_score() const
		{
		if (m_scores_are_evalues)
			return 9999;
		else
			return -9999;
		}
	};
