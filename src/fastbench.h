#pragma once

#include "lookup.h"

#define PARALLEL_SORT	1

class FastBench
	{
public:
	string m_name;
	bool m_scores_are_evalues = false;
	float *m_Scores = 0;
	float m_Sum3 = FLT_MAX;
	float m_SEPQ0_1 = FLT_MAX;
	float m_SEPQ1 = FLT_MAX;
	float m_SEPQ10 = FLT_MAX;
	lookup *m_look = 0;
	uint8_t *m_dope = 0;
	uint32_t m_dope_nhit = 0;
	uint32_t *m_dope_ks = 0;
	uint m_ndom = UINT_MAX;
	uint m_npair = UINT_MAX;

	vector<string> m_Labels;

#if PARALLEL_SORT
	uint *m_ScoreOrder = 0;
	uint  m_ScoreOrderCap = 0;
	bool  m_ScoreOrderInitialized = false;
#else
	uint *m_ScoreOrder = 0;
#endif

public:
	FastBench()
		{
		}

	~FastBench()
		{
		myfree(m_Scores);
		myfree(m_ScoreOrder);
		}

public:
	virtual void SubclassClearHitsAndResults() {}
	virtual void SubclassAppendHit(uint i, uint j, float Score) {}

public:
	void Alloc();
	void ClearHitsAndResults();
	void ReadLookup(const string &FN);
	void ReadDope(const string &FN);
	void SetLookupFromLabels();
	void AppendHit(uint i, uint j, float Score);
	double Bench(const string &Msg = "");
	void SetScoreOrder();
	void SetScoreOrder_Serial();
#if PARALLEL_SORT
	void SetScoreOrder_Parallel();
#endif
	void ReadHits(
		const string &FN,
		uint qidx,
		uint tidx,
		uint scoreidx);
	void WriteHits(
		const string &FN,
		bool IncludeSelf = false,
		bool UpperTriangleOnly = false) const;
	void ReadBits(const string &FN);
	void WriteBits(const string &FN) const;
	bool IsTP(uint LabelIdx_i, uint LabelIdx_j) const;
	bool IsIgnored(uint LabelIdx_i, uint LabelIdx_j) const;
	void log_dope_ks() const;
	bool in_dope(uint k) const
		{
		if (m_dope == 0) return true;
		byte b = m_dope[k/8];
		return b & (1 << k%8);
		}
	bool in_dope(uint i, uint j) const
		{
		if (m_dope == 0) return true;
		uint k = triangle_ij_to_k2(i, j, m_ndom);
		return in_dope(k);
		}
	float get_missing_score() const
		{
		if (m_scores_are_evalues)
			return 9999;
		else
			return -9999;
		}
	};

void guess_fields(
		const string &hitsfn,
		uint &qfi, uint &tfi, uint &sfi,
		bool &scores_are_evalues);
