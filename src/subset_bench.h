#pragma once

#include "scop40bench.h"
#include "xdpmem.h"

typedef void BYTE_SEQ_FN(
	const PDBChain &Chain,
	uint AlphaSize,
	vector<byte> &ByteSeq);

typedef float (*ALIGN_FN)(
	XDPMem &Mem,
	uint LQ, uint LT,
	float Open, float Ext,
	float * const * SWMx);

class SubsetBench
	{
public:
	SBSCORE m_SBS = SBS_Evalue;
	vector<string> m_Labels;
	vector<string> m_Doms;
	vector<string> m_SFs;
	vector<uint> m_SFIdxs;
	vector<uint> m_SFIdxToSize;
	map<string, uint> m_DomToIdx;
	map<string, uint> m_SFToIdx;
	set<uint> m_DopeDomIdxs;
	float m_Open = -999;
	float m_Ext = -999;
	ALIGN_FN m_AF = 0;

//////////////////////////////
// Dope
//////////////////////////////
	uint m_DopeSize = 0;
	uint16_t *m_DomIdxQs = 0;
	uint16_t *m_DomIdxTs = 0;
	bool *m_TPs = 0;
	uint m_NT = 0;

//////////////////////////////
// Hits
//////////////////////////////
	atomic<uint> m_NextPairIdx = 0;
	float *m_Scores = 0;
	uint *m_ScoreOrder = 0;

//////////////////////////////
// Per-feature
//////////////////////////////
	vector<uint> m_AlphaSizes;
	vector<float> m_Weights;
	vector<string> m_ByteSeqFNs;
	vector<vector<vector<byte> > > m_ByteSeqVec;
	vector<string> m_SubstMxFNs;
	vector<float *> m_SubstMxPtrs;

//////////////////////////////
// Per-Dom
//////////////////////////////
	vector<uint> m_Ls;
	uint m_MaxL = 0;

//////////////////////////////
// Per-thread
//////////////////////////////
	uint m_ThreadCount = 0;
	uint m_MaxSeqLength = 0;
	vector<float **> m_SWMxs;
	vector<XDPMem *> m_Mems;

//////////////////////////////
// Bench
//////////////////////////////
	float m_Sum3 = 0;

public:
	uint GetFeatureCount() const { return SIZE(m_AlphaSizes); }
	uint GetDomCount() const { return SIZE(m_Doms); }
	void ReadLookup(const string &FN);
	void AddDom(const string &Dom, const string &ScopId);
	void AllocDope(uint DopeSize);
	void AllocHits();
	void MakeDopeFromHits(const string &FN);
	uint GetDomIdx(const string &Label, bool ErrOk) const;
	uint GetSFIdx(uint DomIdx) const;
	void WriteDope(const string &FN) const;
	void ReadDope(const string &FN);
	void MakeByteSeqs(const string &ChainsFN, BYTE_SEQ_FN BSFn,
		uint AlphaSize, const string &BSFN) const;
	void ReadByteSeqs(const string &FN, uint AS,
		vector<vector<byte> > &ByteSeqs) const;
	void LoadByteSeqs(const vector<string> &FNs);
	void ByteSeqsToFasta(vector<vector<byte> > &ByteSeqs,
		const string &FN) const;
	float * const * MakeSWMx(uint ThreadIdx,
		uint DomIdxQ, uint DomIdxT);
	void ThreadBody(uint ThreadIdx);
	void Search(ALIGN_FN AF);
	float *ReadScoreMx(const string &FN, uint &AlphaSize) const;
	void LoadScoreMxs(vector<string> &FNs);
	void SetWeights(const vector<float> &Weights);
	void Bench(const string &Msg = "");
	void SetScoreOrder();
	float **AllocSWMx() const;

public:
	static void StaticThreadBody(SubsetBench *SB, uint ThreadIdx);
	};

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Weights);
