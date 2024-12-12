#include "myutils.h"
#include "dbsearcher.h"
#include "chainer.h"
#include "mukmerfilter.h"
#include "timing.h"

extern int8_t IntScoreMx_Mu[36][36];

void MuKmerFilter::SetParams(const DSSParams &Params)
	{
	uint k = GetPatternOnes(Params.m_PatternStr);
	asserta(k >= 1 && k < 6);
	m_DictSize = myipow(36, k);
	m_MinHSPScore = Params.m_MKF_MinHSPScore;
	m_X = Params.m_MKF_X;
	}

int MuKmerFilter::MuXDrop(int PosQ, int LQ, int PosT, int LT, int X,
						int &Loi, int &Loj, int &Len)
	{
	StartTimer(MuXDrop);
	Loi = PosQ;
	Loj = PosT;
	Len = 0;
	int i = PosQ;
	int j = PosT;
	const byte *Q = m_ptrMuLettersQ->data();
	const byte *T = m_ptrMuLettersT->data();

	int FwdScore = 0;
	int BestFwdScore = 0;
	int FwdLen = 0;
	while (i < LQ && j < LT)
		{
		byte q = Q[i++];
		byte t = T[j++];
		FwdScore += IntScoreMx_Mu[q][t];
		if (FwdScore > BestFwdScore)
			{
			FwdLen = i - PosQ;
			BestFwdScore = FwdScore;
			}
		else if (FwdScore + X < BestFwdScore)
			break;
		}

	int RevScore = 0;
	int BestRevScore = 0;
	int RevLen = 0;
	i = PosQ - 1;
	j = PosT - 1;
	while (i >= 0 && j >= 0)
		{
		byte q = Q[i];
		byte t = T[j];
		RevScore += IntScoreMx_Mu[q][t];
		if (RevScore > BestRevScore)
			{
			BestRevScore = RevScore;
			Loi = i;
			Loj = j;
			RevLen = PosQ - i;
			}
		else if (RevScore + X < BestRevScore)
			break;
		--i;
		--j;
		}
	int BestScore = BestFwdScore + BestRevScore;
	Len = FwdLen + RevLen;
	EndTimer(MuXDrop);
#if DEBUG
	{
	int CheckScore = 0;
	int i = Loi;
	int j = Loj;
	for (int k = 0; k < Len; ++k)
		{
		asserta(i >= 0 && j >= 0 && i < LQ && j < LT);
		byte q = Q[i++];
		byte t = T[j++];
		CheckScore += IntScoreMx_Mu[q][t];
		}
	asserta(CheckScore = BestScore);
	}
#endif
	return BestScore;
	}

void MuKmerFilter::MuKmerResetQ()
	{
	if (m_KmerHashTableQ == 0)
		{
		m_KmerHashTableQ = myalloc(uint16_t, m_DictSize);
		memset(m_KmerHashTableQ, 0xff, m_DictSize*sizeof(uint16_t));
		}

	if (m_ptrMuKmersQ != 0)
		{
		const uint KmerCountQ = SIZE(*m_ptrMuKmersQ);
		for (uint PosQ = 0; PosQ < KmerCountQ; ++PosQ)
			{
			uint Kmer = (*m_ptrMuKmersQ)[PosQ];
			m_KmerHashTableQ[Kmer] = 0xffff;
			}
		}
#if DEBUG
	{
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		asserta(m_KmerHashTableQ[Kmer] == 0xffff);
	}
#endif
	}

void MuKmerFilter::MuKmerSetQ(const PDBChain &ChainQ,
							  const vector<byte> *ptrMuLettersQ,
							  const vector<uint> *ptrMuKmersQ)
	{
	StartTimer(MuKmerSetQ);
	m_ChainQ = &ChainQ;

	m_ptrMuLettersQ = ptrMuLettersQ;
	m_ptrMuKmersQ = ptrMuKmersQ;

	const uint KmerCountQ = SIZE(*m_ptrMuKmersQ);
	asserta(KmerCountQ < 0xffff);
	for (uint PosQ = 0; PosQ < KmerCountQ; ++PosQ)
		{
		uint Kmer = (*m_ptrMuKmersQ)[PosQ];
		assert(Kmer < m_DictSize);
		if (m_KmerHashTableQ[Kmer] == 0xffff)
			m_KmerHashTableQ[Kmer] = PosQ;
		}
	EndTimer(MuKmerSetQ);
	}

bool MuKmerFilter::MuKmerAln(const PDBChain &ChainT,
						   const vector<byte> &MuLettersT,
						   const vector<uint> &MuKmersT)
	{
	StartTimer(MuKmerAln);
	m_ChainT = &ChainT;
	m_ptrMuLettersT = &MuLettersT;
	++m_MuKmerFilterPairCount;
	const uint KmerCountT = SIZE(MuKmersT);
	int LQ = int(m_ChainQ->GetSeqLength());
	int LT = int(m_ChainT->GetSeqLength());
	bool FoundHSP = false;
	int BestHSPScore = 0;
	m_MuKmerHSPLois.clear();
	m_MuKmerHSPLojs.clear();
	m_MuKmerHSPLens.clear();
	m_MuKmerHSPScores.clear();
	m_BestChainScore = 0;
	for (uint PosT = 0; PosT < KmerCountT; ++PosT)
		{
		uint KmerT = MuKmersT[PosT];
		assert(KmerT < m_DictSize);
		uint PosQ = m_KmerHashTableQ[KmerT];
		if (PosQ != 0xffff)
			{
			EndTimer(MuKmerAln);
			int Loi, Loj, Len;
			int Score = MuXDrop(int(PosQ), LQ, int(PosT), LT, m_X, Loi, Loj, Len);
			StartTimer(MuKmerAln);
			if (Score >= m_MinHSPScore)
				{
				if (Score > BestHSPScore)
					{
					FoundHSP = true;
					bool OldHSP = false;
					for (uint i = 0; i < SIZE(m_MuKmerHSPLois); ++i)
						{
						if (Loi == m_MuKmerHSPLois[i])
							{
							OldHSP = true;
							break;
							}
						}
					if (!OldHSP)
						{
						m_MuKmerHSPLois.push_back(Loi);
						m_MuKmerHSPLojs.push_back(Loj);
						m_MuKmerHSPLens.push_back(Len);
						m_MuKmerHSPScores.push_back(Score);
						}
					}
				}
			}
		}
	EndTimer(MuKmerAln);
	if (FoundHSP)
		{
		++m_MuKmerFilterHitCount;
		ChainHSPs();
		}
	return FoundHSP;
	}

void MuKmerFilter::ChainHSPs()
	{
	const uint N = SIZE(m_MuKmerHSPLois);
	asserta(SIZE(m_MuKmerHSPLens) == N);
	asserta(SIZE(m_MuKmerHSPScores) == N);
	vector<uint> Los;
	vector<uint> His;
	vector<float> Scores;
	for (uint i = 0; i < N; ++i)
		{
		uint Lo = m_MuKmerHSPLois[i];
		uint Hi = Lo + m_MuKmerHSPLens[i] - 1;
		Los.push_back(Lo);
		His.push_back(Hi);
		Scores.push_back((float) m_MuKmerHSPScores[i]);
		}
	Chainer C;
	vector<uint> Idxs;
	m_BestChainScore = (int) C.Chain(Los, His, Scores, Idxs);
	m_ChainLo_i = 0;
	m_ChainHi_i = 0;
	m_ChainLo_j = 0;
	m_ChainHi_j = 0;
	const uint M = SIZE(Idxs);
	if (M == 0)
		return;
	for (uint k = 0; k < M; ++k)
		{
		int Idx = Idxs[k];
		int Loi = m_MuKmerHSPLois[Idx];
		int Loj = m_MuKmerHSPLojs[Idx];
		int Hii = Loi + m_MuKmerHSPLens[Idx] - 1;
		int Hij = Loj + m_MuKmerHSPLens[Idx] - 1;
		if (k == 0)
			{
			m_ChainLo_i = Loi;
			m_ChainHi_i = Hii;

			m_ChainLo_j = Loj;
			m_ChainHi_j = Hij;
			}
		else
			{
			m_ChainLo_i = min(Loi, m_ChainLo_i);
			m_ChainHi_i = max(Hii, m_ChainHi_i);

			m_ChainLo_j = min(Loj, m_ChainLo_j);
			m_ChainHi_j = max(Hij, m_ChainHi_j);
			}
		}
	}
