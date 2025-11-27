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
	vector<Paralign *> m_PAs;
	vector<uint> m_QueryIdxs;

	SCOP40Bench m_SB;

	string m_AlignMethod;
	string m_SubstMxName;
	string m_ByteSeqMethod;

	uint m_SeqCount = UINT_MAX;
	uint m_PairCount = UINT_MAX;

public:
	void GetByteSeqs_numu(const string &FN);
	void GetByteSeqs_muletters(const string &FN);
	void GetByteSeqs_dss3(const string &FN);
	void GetByteSeqs(const string &FN, const string &Method);
	void Search(const string &AlignMethod, string SubstMxName);
	void SetQuery(uint ThreadIdx, uint i);
	void Align(uint ThreadIdx, uint i, uint j);
	void AppendHit(uint i, uint j, float Score);
	void Bench(const string &LookupFN);
	void WriteHits(const string &FN) const;
	void SetDomIdxs();
	void ReadLookup(const string &FN)
		{
		m_SB.ReadLookup(FN);
		}
	};

void FixMuByteSeq(vector<byte> &ByteSeq);