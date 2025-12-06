#pragma once

#include "paralign.h"

class PDBChain;

class ParaSearch
	{
public:
	vector<PDBChain *> m_Chains;
	vector<string> m_Labels;
	vector<vector<byte> > m_ByteSeqs;
	vector<uint> m_DomIdxs;

	float *m_Scores = 0;
	float m_Sum3 = FLT_MAX;

	vector<Paralign *> m_PAs;
	vector<uint> m_QueryIdxs;

	string m_AlignMethod;
	string m_SubstMxName;
	string m_ByteSeqMethod;

	uint m_SeqCount = UINT_MAX;
	uint m_PairCount = UINT_MAX;
	vector<uint> m_LabelIdxToSFIdx;
	vector<string> m_SFs;
	vector<uint> m_SFIdxToSize;
	unordered_map<string, uint> m_SFToIdx;
	uint m_NT;	// upper triangle only
	uint m_NF;	// upper triangle only
	uint *m_ScoreOrder;

public:
	static vector<FEATURE> m_NuFs;

public:
	void ClearHitsAndResults();
	void SetGapParams(int Open, int Ext);
	void GetByteSeqs(const string &FN, const string &Method);
	void Search(const string &AlignMethod);
	void SetQuery(uint ThreadIdx, uint i);
	void Align(uint ThreadIdx, uint i, uint j);
	void AppendHit(uint i, uint j, float Score);
	void Bench();
	void SetScoreOrder();
	void WriteHits(const string &FN) const;
	void SetLookupFromLabels();
	void AddDom(const string &Dom, const string &SF, uint LabelIdx);
	void MakeSubset(ParaSearch &Subset, uint SubsetPct);

	private:
	void GetByteSeqs_nu(const string &FN);
	void GetByteSeqs_numu(const string &FN);
	void GetByteSeqs_muletters(const string &FN);
	void GetByteSeqs_dss3(const string &FN);

	};

void FixMuByteSeq(vector<byte> &ByteSeq);