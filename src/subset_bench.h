#pragma once

#include "scop40bench.h"

typedef void BYTE_SEQ_FN(const PDBChain &Chain,
	vector<byte> &ByteSeq);

class SubsetBench
	{
public:
	SBSCORE m_SBS = SBS_Evalue;
	vector<string> m_Labels;
	vector<string> m_Doms;
	vector<string> m_SFs;
	vector<uint> m_SFIdxs;
	map<string, uint> m_DomToIdx;
	map<string, uint> m_SFToIdx;
	set<uint> m_DopeDomIdxs;

//////////////////////////////
// Dope
//////////////////////////////
	uint m_DopeSize = 0;
	uint16_t *m_DomIdxQs = 0;
	uint16_t *m_DomIdxTs = 0;
	bool *m_TPs = 0;
//////////////////////////////

public:
	void ReadLookup(const string &FN);
	void AddDom(const string &Dom, const string &ScopId);
	void AllocDope(uint DopeSize);
	void MakeDopeFromHits(const string &FN);
	uint GetDomIdx(const string &Label) const;
	uint GetSFIdx(uint DomIdx) const;
	void WriteDope(const string &FN) const;
	void ReadDope(const string &FN);
	void MakeByteSeqs(const string &ChainsFN,
		BYTE_SEQ_FN BSFn, const string &BSFN);
	};