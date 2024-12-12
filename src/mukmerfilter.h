#pragma once

#include "dss.h"

class MuKmerFilter
	{
private:
	const DSSParams *m_Params = 0;

public:
	uint m_DictSize = 0;
	int m_MinHSPScore = INT_MAX;
	int m_X = INT_MAX;
	const vector<byte> *m_ptrMuLettersQ = 0;
	const vector<byte> *m_ptrMuLettersT = 0;
	const vector<uint> *m_ptrMuKmersQ = 0;
	uint16_t *m_KmerHashTableQ = 0;
	uint m_MuKmerFilterPairCount = 0;
	uint m_MuKmerFilterHitCount = 0;
	const PDBChain *m_ChainQ = 0;
	const PDBChain *m_ChainT = 0;
	vector<int> m_MuKmerHSPLois;
	vector<int> m_MuKmerHSPLojs;
	vector<int> m_MuKmerHSPLens;
	vector<int> m_MuKmerHSPScores;
	int m_BestChainScore = 0;
	int m_ChainLo_i = 0;
	int m_ChainHi_i = 0;
	int m_ChainLo_j = 0;
	int m_ChainHi_j = 0;

public:
	void SetParams(const DSSParams &Params);
	void MuKmerSetQ(const PDBChain &ChainQ,
					const vector<byte> *ptrMuLettersQ,
					const vector<uint> *ptrMuKmersQ);
	void MuKmerResetQ();
	bool MuKmerAln(const PDBChain &ChainT,
				   const vector<byte> &MuLettersT,
				   const vector<uint> &MuKmersT);
	int MuXDrop(int PosQ, int LQ, int PosT, int LT, int X,
				int &Loi, int &Loj, int &Len);
	void ChainHSPs();
	};
