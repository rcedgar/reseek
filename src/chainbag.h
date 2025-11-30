#pragma once

class PDBChain;

class ChainBag
	{
public:
	const PDBChain *m_ptrChain = 0;
	const vector<vector<byte> > *m_ptrProfile = 0;
	const vector<byte> *m_ptrMuLetters = 0;
	const vector<uint> *m_ptrMuKmers = 0;
	const void *m_ptrProfPara8 = 0;
	const void *m_ptrProfPara16 = 0;
	const void *m_ptrProfParaRev8 = 0;
	const void *m_ptrProfParaRev16 = 0;
	const uint16_t *m_ptrKmerHashTableQ = 0;
	float m_SelfRevScore = FLT_MAX;

public:
	void Validate(const char *Msg) const;
	};
