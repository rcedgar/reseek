#include "myutils.h"
#include "dbsearcher.h"
#include "timing.h"

//static const uint DictSize = 36*36*36*36*36;
//static const string PatternStr = "11111";
//static const uint DictSize = 36*36*36*36;
//static const string PatternStr = "1111";
static const uint DictSize = 36*36*36;
static const string PatternStr = "111";
extern int8_t IntScoreMx_Mu[36][36];

#if MUKMERS
int DBSearcher::MuXDrop(int PosQ, int LQ, int PosT, int LT, int X)
	{
	StartTimer(MuXDrop);
	int i = PosQ;
	int j = PosT;
	const byte *Q = m_MuLettersQ.data();
	const byte *T = (*m_ptrMuLettersT).data();

	int FwdScore = 0;
	int BestFwdScore = 0;
	while (i < LQ && j < LT)
		{
		byte q = Q[i++];
		byte t = T[j++];
		FwdScore += IntScoreMx_Mu[q][t];
		if (FwdScore > BestFwdScore)
			BestFwdScore = FwdScore;
		else if (FwdScore + X < BestFwdScore)
			break;
		}

	int RevScore = 0;
	int BestRevScore = 0;
	i = PosQ - 1;
	j = PosT - 1;
	while (i >= 0 && j >= 0)
		{
		byte q = Q[i--];
		byte t = T[j--];
		RevScore += IntScoreMx_Mu[q][t];
		if (RevScore > BestRevScore)
			BestRevScore = RevScore;
		else if (RevScore + X < BestRevScore)
			break;
		}
	int BestScore = BestFwdScore + BestRevScore;
	EndTimer(MuXDrop);
	return BestScore;
	}

void DBSearcher::MuKmerResetQ()
	{
	if (m_KmerHashTableQ == 0)
		{
		asserta(m_Params->m_PatternStr == PatternStr);
		m_KmerHashTableQ = myalloc(uint16_t, DictSize);
		memset(m_KmerHashTableQ, 0xff, DictSize*sizeof(uint16_t));
		}

	const uint KmerCountQ = SIZE(m_MuKmersQ);
	for (uint PosQ = 0; PosQ < KmerCountQ; ++PosQ)
		{
		uint Kmer = m_MuKmersQ[PosQ];
		m_KmerHashTableQ[Kmer] = 0xffff;
		}
#if DEBUG
	{
	for (uint Kmer = 0; Kmer < DictSize; ++Kmer)
		asserta(m_KmerHashTableQ[Kmer] == 0xffff);
	}
#endif
	}

void DBSearcher::MuKmerSetQ(const PDBChain &ChainQ)
	{
	StartTimer(MuKmerSetQ);
	m_ChainQ = &ChainQ;

	m_D.Init(ChainQ);
	m_D.GetMuLetters(m_MuLettersQ);
	m_D.GetMuKmers(m_MuLettersQ, m_MuKmersQ);

	const uint KmerCountQ = SIZE(m_MuKmersQ);
	asserta(KmerCountQ < 0xffff);
	for (uint PosQ = 0; PosQ < KmerCountQ; ++PosQ)
		{
		uint Kmer = m_MuKmersQ[PosQ];
		assert(Kmer < DictSize);
		if (m_KmerHashTableQ[Kmer] == 0xffff)
			m_KmerHashTableQ[Kmer] = PosQ;
		}
	EndTimer(MuKmerSetQ);
	}

bool DBSearcher::MuKmerAln(const PDBChain &ChainT, double Evalue,
						   const vector<byte> &MuLettersT,
						   const vector<uint> &MuKmersT)
	{
	StartTimer(MuKmerAln);
	m_ChainT = &ChainT;
	m_ptrMuLettersT = &MuLettersT;
	++m_QTPairCount;
	if (Evalue <= 1)
		++m_MuEle1Count;
	if (Evalue <= 10)
		++m_MuEle10Count;
	const uint KmerCountT = SIZE(MuKmersT);
	int LQ = int(m_ChainQ->GetSeqLength());
	int LT = int(m_ChainT->GetSeqLength());
	bool FoundHSP = false;
	for (uint PosT = 0; PosT < KmerCountT; ++PosT)
		{
		uint KmerT = MuKmersT[PosT];
		assert(KmerT < DictSize);
		uint PosQ = m_KmerHashTableQ[KmerT];
		if (PosQ != 0xffff)
			{
			EndTimer(MuKmerAln);
			++m_QTSeedCount;
			int Score = MuXDrop(int(PosQ), LQ, int(PosT), LT, 10);
			StartTimer(MuKmerAln);
			if (Score >= 50)
				{
				FoundHSP = true;
				++m_MuHSPCount;
				if (Evalue <= 1)
					++m_MuHSPCount1;
				if (Evalue <= 10)
					++m_MuHSPCount10;
				else
					++m_MuFPHSPCount;
				break;
				}
			}
		}
	EndTimer(MuKmerAln);
	return FoundHSP;
	}
#endif // MUKMERS