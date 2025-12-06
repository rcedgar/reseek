#pragma once

#include "paralign.h"
#include "scop40bench.h"

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

//	SCOP40Bench m_SB;

	string m_AlignMethod;
	string m_SubstMxName;
	string m_ByteSeqMethod;

	uint m_SeqCount = UINT_MAX;
	uint m_PairCount = UINT_MAX;

public:
	static vector<FEATURE> m_NuFs;
	static vector<uint> m_LabelIdxToSFIdx;
	static vector<string> m_SFs;
	static vector<uint> m_SFIdxToSize;
	static unordered_map<string, uint> m_SFToIdx;
	static uint m_NT;	// upper triangle only
	static uint m_NF;	// upper triangle only
	static uint *m_ScoreOrder;

public:
	void ClearHitsAndResults();
	void SetGapParams(int Open, int Ext);
	void GetByteSeqs_nu(const string &FN);
	void GetByteSeqs_numu(const string &FN);
	void GetByteSeqs_muletters(const string &FN);
	void GetByteSeqs_dss3(const string &FN);
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
	};

void FixMuByteSeq(vector<byte> &ByteSeq);