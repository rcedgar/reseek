#pragma once

#include "fastbench.h"
#include "paralign.h"

class PDBChain;

class ParaSearch : public FastBench
	{
public:
	vector<PDBChain *> m_Chains;
	vector<vector<byte> > m_ByteSeqs;
	vector<uint> m_DomIdxs;

	float *m_SelfScores_rev = 0;
	float *m_Scores_fwd = 0;
	float *m_Scores_rev = 0;

	vector<Paralign *> m_PAs;
	vector<uint> m_QueryIdxs;

	string m_AlignMethod;
	string m_SubstMxName;
	string m_ByteSeqMethod;

	bool m_DoReverse = false;

public:
	static vector<FEATURE> m_NuFs;

public:
	virtual void SubclassClearHitsAndResults();
	virtual void SubclassAppendHit(uint i, uint j, float Score);

public:
	void ClearHitsAndResults();
	void SetGapParams(int Open, int Ext);
	void GetByteSeqs(const string &FN, const string &Method);
	void Search(const string &AlignMethod, bool DoReverse);
	void InitThreads(const string &AlignMethod, bool DoReverse);
	void SetQuery(uint ThreadIdx, uint i);
	void Align(uint ThreadIdx, uint i, uint j);
	void WriteRevTsv(const string &FN) const;
	void MakeSubset(ParaSearch &Subset, uint SubsetPct);

	void SetSelfScores_rev(const string &AlignMethod);
	float GetSelfScore_rev(Paralign &PA, uint ChainIdx);
	void AppendHit_rev(uint i, uint j, float Score);
	void BenchRev(const string &Msg, float SelfWeight, float RevWeight);

private:
	void GetByteSeqs_nu(const string &FN);
	void GetByteSeqs_numu(const string &FN);
	void GetByteSeqs_muletters(const string &FN);
	void GetByteSeqs_dss3(const string &FN);
	};

void FixMuByteSeq(vector<byte> &ByteSeq);