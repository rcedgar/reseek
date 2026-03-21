#pragma once

#include "lookup.h"

class FastBench
	{
public:
	float *m_Scores = 0;
	float m_Sum3 = FLT_MAX;
	uint m_SeqCount = UINT_MAX;
	uint m_PairCount = UINT_MAX;
	uint *m_ScoreOrder = 0;
	vector<string> m_Labels;
	lookup *m_look = 0;

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
	void SetLookupFromLabels();
	void AppendHit(uint i, uint j, float Score);
	void Bench(const string &Msg = "");
	void SetScoreOrder();
	void WriteHits(const string &FN, bool IncludeSelf = false) const;
	bool IsTP(uint LabelIdx_i, uint LabelIdx_j) const;
	};
