#pragma once

#include "dss.h"
#include "chainer.h"

const uint HASHW = 4;

class MuKmerFilter
	{
private:
	const DSSParams *m_Params = 0;

public:
	uint m_DictSize = 0;
	const vector<byte> *m_ptrMuLettersQ = 0;
	const vector<byte> *m_ptrMuLettersT = 0;
	const vector<uint> *m_ptrMuKmersQ = 0;
	uint16_t *m_KmerHashTableQ = 0;
	vector<int> m_MuKmerHSPLois;
	vector<int> m_MuKmerHSPLojs;
	vector<int> m_MuKmerHSPLens;
	vector<int> m_MuKmerHSPScores;
	int m_BestChainScore = 0;
	int m_ChainLo_i = 0;
	int m_ChainHi_i = 0;
	int m_ChainLo_j = 0;
	int m_ChainHi_j = 0;
	vector<int> m_ChainHSPLois;
	vector<int> m_ChainHSPLojs;
	vector<int> m_ChainHSPLens;
	vector<float> Scores;

public:
	Chainer m_C;

public:
	~MuKmerFilter()
		{
		myfree(m_KmerHashTableQ);
		}

public:
	static uint m_PairCount;
	static uint m_PairWithHSPCount;
	static uint m_BoxAlnCount;
	static uint m_ChainCount;
	static uint m_ChainHSPCount;
	static double m_SumBoxFract;
	static mutex m_Lock;

public:
	void SetParams(const DSSParams &Params);
	void SetQ(const vector<byte> *ptrMuLettersQ,
					const vector<uint> *ptrMuKmersQ);

// TODO -- maybe no need to support ResetQ()?
	void ResetQ();
	void Align(const vector<byte> &MuLettersT,
				   const vector<uint> &MuKmersT);
	int MuXDrop(int PosQ, int LQ, int PosT, int LT, int X,
				int &Loi, int &Loj, int &Len);
	void ChainHSPs();
	uint GetQL() const { return SIZE(*m_ptrMuLettersQ); };
	uint GetTL() const { return SIZE(*m_ptrMuLettersT); };

public:
	static void Stats();
	};
