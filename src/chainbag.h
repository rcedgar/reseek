#pragma once

class PDBChain;

class ChainBag
	{
public:
	const PDBChain *m_ptrChain = 0;
	const vector<vector<byte> > *m_ptrProfile = 0;
	const vector<byte> *m_ptrMuLetters = 0;
	const vector<uint> *m_ptrMuKmers = 0;
	void *m_ptrProfPara = 0;
	void *m_ptrProfParaRev = 0;
	uint16_t *m_ptrKmerHashTableQ = 0;
	float m_ptrSelfRevScore = FLT_MAX;
	};
