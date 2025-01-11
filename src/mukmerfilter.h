#pragma once

#include "dss.h"
#include "chainer.h"
#include "chainbag.h"

const uint HASHW = 4;

class MuKmerFilter
	{
private:
	const DSSParams *m_Params = 0;

	const ChainBag *m_ptrBagQ = 0;
	const vector<byte> *m_ptrMuLettersQ = 0;
	const vector<uint> *m_ptrMuKmersQ = 0;
	uint16_t *m_ptrKmerHashTableQ = 0;

	const vector<byte> *m_ptrMuLettersT = 0;

public:
	string m_LabelQ;
	uint m_DictSize = 0;
	vector<int> m_MuKmerHSPLois;
	vector<int> m_MuKmerHSPLojs;
	vector<int> m_MuKmerHSPLens;
	vector<int> m_MuKmerHSPScores;
	int m_BestChainScore = 0;
	int m_BestHSPScore = 0;
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
		myfree(m_ptrKmerHashTableQ);
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
	const uint16_t *GetHashTableQ() const
		{
		return m_ptrKmerHashTableQ;
		}
	void SetParams(const DSSParams &Params);
	void ResetQ();
	void SetQ(const string &LabelQ,
			  const vector<byte> *ptrMuLettersQ,
			  const vector<uint> *ptrMuKmersQ);
	void SetBagQ(const ChainBag &BagQ);
	int GetMaxHSPScore(const vector<byte> &MuLettersT,
					   const vector<uint> &MuKmersT);
	void Align(const vector<byte> &MuLettersT,
				   const vector<uint> &MuKmersT);
	void AlignBag(const ChainBag &BagT);
	int MuXDrop(int PosQ, int LQ, int PosT, int LT, int X,
				int &Loi, int &Loj, int &Len) const;
	void ChainHSPs();
	uint GetQL() const { return SIZE(*m_ptrMuLettersQ); };
	uint GetTL() const { return SIZE(*m_ptrMuLettersT); };
#if DEBUG
	void Validate() const;
	void Dump(FILE *f, const char *Msg) const;
#endif

public:
	static void Stats();
	};
