#pragma once

class FastBench
	{
public:
	float *m_Scores = 0;
	float m_Sum3 = FLT_MAX;
	uint m_SeqCount = UINT_MAX;
	uint m_PairCount = UINT_MAX;
	vector<string> m_Labels;
	vector<uint> m_LabelIdxToSFIdx;
	vector<string> m_SFs;
	vector<uint> m_SFIdxToSize;
	unordered_map<string, uint> m_SFToIdx;
	uint m_NT = UINT_MAX;
	uint m_NF = UINT_MAX;
	uint *m_ScoreOrder = 0;

public:
	FastBench()
		{
		m_Labels.reserve(12000);
		m_LabelIdxToSFIdx.reserve(12000);
		m_LabelIdxToSFIdx.reserve(12000);
		m_SFs.reserve(2000);
		m_SFIdxToSize.reserve(2000);
		m_SFToIdx.reserve(2000);
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
	void AppendLabel(const string &Label);
	void Alloc();
	void ClearHitsAndResults();
	void AddDom(const string &Dom, const string &SF, uint LabelIdx);
	void AppendHit(uint i, uint j, float Score);
	void Bench(const string &Msg = "");
	void SetScoreOrder();
	void WriteHits(const string &FN) const;
	void SetLookupFromLabels();
	};
