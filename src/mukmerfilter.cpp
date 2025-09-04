#include "myutils.h"
#include "dbsearcher.h"
#include "chainer.h"
#include "mukmerfilter.h"
#include "mudex.h"
#include "alpha.h"
#include "timing.h"

extern int8_t IntScoreMx_Mu[36][36];

uint MuKmerFilter::m_PairCount;
uint MuKmerFilter::m_PairWithHSPCount;
uint MuKmerFilter::m_BoxAlnCount;
uint MuKmerFilter::m_ChainCount;
uint MuKmerFilter::m_ChainHSPCount;
double MuKmerFilter::m_SumBoxFract;
mutex MuKmerFilter::m_Lock;

#if DEBUG
void MuKmerFilter::Dump(FILE *f, const char *Msg) const
	{
	if (f == 0)
		return;
	fprintf(f, "MuKmerFilter::Dump(%s) m_Label=%s\n", Msg, m_LabelQ.c_str());
	fprintf(f, "Letters %p ", m_ptrMuLettersQ);
	if (m_ptrMuLettersQ != 0)
		{
		fprintf(f, "[%u]", SIZE(*m_ptrMuLettersQ));
		for (uint i = 0; i < SIZE(*m_ptrMuLettersQ); ++i)
			fprintf(f, " %u", (*m_ptrMuLettersQ)[i]);
		}
	fprintf(f, "\n");

	fprintf(f, "Kmers %p ", m_ptrMuKmersQ);
	if (m_ptrMuKmersQ != 0)
		{
		const uint n = SIZE(*m_ptrMuKmersQ);
		fprintf(f, "[%u]", n);
		for (uint i = 0; i < n; ++i)
			fprintf(f, " %u", (*m_ptrMuKmersQ)[i]);
		}
	fprintf(f, "\n");

	uint Busy = 0;
	for (uint Kmer = 0; Kmer < m_DictSize*HASHW; ++Kmer)
		{
		if (m_ptrKmerHashTableQ[Kmer] != 0xffff)
			++Busy;
		}
	Log("Busy=%u\n", Busy);
	}

void MuKmerFilter::Validate() const
	{
	asserta(m_ptrKmerHashTableQ != 0);
	if (m_ptrMuLettersQ == 0)
		{
		asserta(m_ptrMuKmersQ == 0);
		for (uint Kmer = 0; Kmer < m_DictSize*HASHW; ++Kmer)
			asserta(m_ptrKmerHashTableQ[Kmer] == 0xffff);
		return;
		}
	asserta(m_ptrMuKmersQ != 0);
	const uint KmerCountQ = SIZE(*m_ptrMuKmersQ);
	set<uint> KmerSetQ;
	for (uint PosQ = 0; PosQ < KmerCountQ; ++PosQ)
		{
		uint Kmer = (*m_ptrMuKmersQ)[PosQ];
		KmerSetQ.insert(Kmer);
		assert(Kmer < m_DictSize);
		asserta(m_ptrKmerHashTableQ[HASHW*Kmer] != 0xffff);
		uint low = UINT_MAX;
		for (uint w = 1; w < HASHW; ++w)
			{
			if (m_ptrKmerHashTableQ[HASHW*Kmer+w] == 0xffff)
				{
				low = w;
				break;
				}
			}
		if (low != UINT_MAX)
			{
			for (uint w = low+1; w < HASHW; ++w)
				asserta(m_ptrKmerHashTableQ[HASHW*Kmer+w] == 0xffff);
			}
		}
	for (uint Kmer = 0; Kmer < m_DictSize; ++Kmer)
		{
		if (KmerSetQ.find(Kmer) != KmerSetQ.end())
			continue;
		for (uint w = 0; w < HASHW; ++w)
			asserta(m_ptrKmerHashTableQ[HASHW*Kmer+w] == 0xffff);
		}
	}
#endif // DEBUG

void MuKmerFilter::SetParams(const DSSParams &Params)
	{
	uint k = GetPatternOnes(Params.m_MKFPatternStr);
	asserta(k >= 1 && k < 6);
	m_DictSize = myipow(36, k);
	m_Params = &Params;
	}

int MuKmerFilter::MuXDrop(int PosQ, int LQ, int PosT, int LT, int X,
						int &Loi, int &Loj, int &Len) const
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
	asserta(CheckScore == BestScore);
	}
#endif
	return BestScore;
	}

void MuKmerFilter::ResetQ()
	{
	if (m_ptrMuKmersQ != 0)
		{
		asserta(m_ptrMuLettersQ != 0);
		const uint KmerCountQ = SIZE(*m_ptrMuKmersQ);
		for (uint PosQ = 0; PosQ < KmerCountQ; ++PosQ)
			{
			uint Kmer = (*m_ptrMuKmersQ)[PosQ];
			assert(Kmer < m_DictSize);
			for (uint w = 0; w < HASHW; ++w)
				m_ptrKmerHashTableQ[HASHW*Kmer+w] = 0xffff;
			}
		m_ptrMuKmersQ = 0;
		m_ptrMuLettersQ = 0;
		m_LabelQ.clear();
		}
#if DEBUG
	{
	Validate();
	for (uint Kmer = 0; Kmer < m_DictSize*HASHW; ++Kmer)
		asserta(m_ptrKmerHashTableQ[Kmer] == 0xffff);
	}
#endif
	}

void MuKmerFilter::InitHashTable(uint16_t *HT) const
	{
	memset(HT, 0xff, m_DictSize*HASHW*sizeof(uint16_t));
	}

void MuKmerFilter::SetHashTable(const vector<uint> &Kmers, uint16_t *HT) const
	{
	const uint KmerCount = SIZE(Kmers);
	asserta(KmerCount < 0xffff);
	for (uint PosQ = 0; PosQ < KmerCount; ++PosQ)
		{
		uint Kmer = Kmers[PosQ];
		assert(Kmer < m_DictSize);
		for (uint w = 0; w < HASHW; ++w)
			{
			if (HT[Kmer*HASHW+w] == 0xffff)
				{
				HT[Kmer*HASHW+w] = PosQ;
				break;
				}
			}
		}
	}

uint16_t *MuKmerFilter::CreateEmptyHashTable() const\
	{
	uint16_t *HT = myalloc(uint16_t, m_DictSize*HASHW);
	memset(HT, 0xff, m_DictSize*HASHW*sizeof(uint16_t));
	return HT;
	}

void MuKmerFilter::SetQ(const string &LabelQ,
	const vector<byte> *ptrMuLettersQ,
	const vector<uint> *ptrMuKmersQ)
	{
	if (m_ptrKmerHashTableQ == 0)
		{
		m_ptrKmerHashTableQ = CreateEmptyHashTable();
		InitHashTable(m_ptrKmerHashTableQ);
		}
	else
		ResetQ();

#if DEBUG
	Validate();
#endif
	m_LabelQ = LabelQ;
	StartTimer(MuKmerSetQ);
	//m_ChainQ = &ChainQ;

	m_ptrMuLettersQ = ptrMuLettersQ;
	m_ptrMuKmersQ = ptrMuKmersQ;

	//const uint KmerCountQ = SIZE(*m_ptrMuKmersQ);
	//asserta(KmerCountQ < 0xffff);
	//for (uint PosQ = 0; PosQ < KmerCountQ; ++PosQ)
	//	{
	//	uint Kmer = (*m_ptrMuKmersQ)[PosQ];
	//	assert(Kmer < m_DictSize);
	//	for (uint w = 0; w < HASHW; ++w)
	//		{
	//		if (m_ptrKmerHashTableQ[Kmer*HASHW+w] == 0xffff)
	//			{
	//			m_ptrKmerHashTableQ[Kmer*HASHW+w] = PosQ;
	//			break;
	//			}
	//		}
	//	}
	SetHashTable(*m_ptrMuKmersQ, m_ptrKmerHashTableQ);
#if DEBUG
	Validate();
#endif
	EndTimer(MuKmerSetQ);
	}

int MuKmerFilter::GetMaxHSPScore(const vector<byte> &MuLettersT,
							 const vector<uint> &MuKmersT)
	{
	StartTimer(GetMaxHSPScore);
	m_ptrMuLettersT = &MuLettersT;
	const uint KmerCountT = SIZE(MuKmersT);
	int LQ = int(GetQL());
	int LT = int(SIZE(MuLettersT));
	const int X1 = m_Params->m_MKF_X1;
	int MaxHSPScore = 0;
	for (uint PosT = 0; PosT < KmerCountT; ++PosT)
		{
		uint KmerT = MuKmersT[PosT];
		assert(KmerT < m_DictSize);
		for (uint w = 0; w < HASHW; ++w)
			{
			uint PosQ = m_ptrKmerHashTableQ[HASHW*KmerT+w];
			if (PosQ != 0xffff)
				{
				int Loi, Loj, Len;
				int Score = MuXDrop(int(PosQ), LQ, int(PosT), LT, X1, Loi, Loj, Len);
				if (Score > MaxHSPScore)
					MaxHSPScore = Score;
				}
			}
		}
	EndTimer(GetMaxHSPScore);
	return MaxHSPScore;
	}

void MuKmerFilter::AlignBag(const ChainBag &BagT)
	{
	m_ptrBagT = &BagT;
	asserta(BagT.m_ptrMuLetters != 0);
	asserta(BagT.m_ptrMuKmers != 0);
	Align(*BagT.m_ptrMuLetters, *BagT.m_ptrMuKmers);
	}

void MuKmerFilter::Align(const vector<byte> &MuLettersT,
							 const vector<uint> &MuKmersT)
	{
	StartTimer(MuKmerAln);
	m_C.Clear();
	//m_ChainT = &ChainT;
	m_ptrMuLettersT = &MuLettersT;
	m_Lock.lock();
	++m_PairCount;
	m_Lock.unlock();
	const uint KmerCountT = SIZE(MuKmersT);
	int LQ = int(GetQL());
	int LT = int(GetTL());
	m_MuKmerHSPLois.clear();
	m_MuKmerHSPLojs.clear();
	m_MuKmerHSPLens.clear();
	m_MuKmerHSPScores.clear();
	m_ChainHSPLois.clear();
	m_ChainHSPLojs.clear();
	m_ChainHSPLens.clear();
	m_BestChainScore = 0;
	m_BestHSPScore = 0;
	bool FoundHSP = false;
	const int MinHSPScore = m_Params->m_MKF_MinHSPScore;
	const int X1 = m_Params->m_MKF_X1;
	for (uint PosT = 0; PosT < KmerCountT; ++PosT)
		{
		uint KmerT = MuKmersT[PosT];
		assert(KmerT < m_DictSize);
		for (uint w = 0; w < HASHW; ++w)
			{
			uint PosQ = m_ptrKmerHashTableQ[HASHW*KmerT+w];
			if (PosQ != 0xffff)
				{
				EndTimer(MuKmerAln);
				int Loi, Loj, Len;
				int Score = MuXDrop(int(PosQ), LQ, int(PosT), LT, X1, Loi, Loj, Len);
				StartTimer(MuKmerAln);
				if (Score >= MinHSPScore)
					{
					FoundHSP = true;
					if (Score > m_BestHSPScore)
						{
						m_BestHSPScore = Score;
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
		}
	EndTimer(MuKmerAln);
	if (FoundHSP)
		{
		m_Lock.lock();
		++m_PairWithHSPCount;
		m_Lock.unlock();
		ChainHSPs();
		}
	}

void MuKmerFilter::ChainHSPs()
	{
	m_ChainHSPLois.clear();
	m_ChainHSPLojs.clear();
	m_ChainHSPLens.clear();

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
	vector<uint> Idxs;
	m_BestChainScore = (int) m_C.Chain(Los, His, Scores, Idxs);
	m_ChainLo_i = 0;
	m_ChainHi_i = 0;
	m_ChainLo_j = 0;
	m_ChainHi_j = 0;
	const uint M = SIZE(Idxs);
	if (M == 0)
		return;

	m_Lock.lock();
	++m_ChainCount;
	m_ChainHSPCount += M;
	m_Lock.unlock();

	for (uint k = 0; k < M; ++k)
		{
		int Idx = Idxs[k];
		int Loi = m_MuKmerHSPLois[Idx];
		int Loj = m_MuKmerHSPLojs[Idx];
		int Len = m_MuKmerHSPLens[Idx];
		int Hii = Loi + Len - 1;
		int Hij = Loj + Len - 1;

		m_ChainHSPLois.push_back(Loi);
		m_ChainHSPLojs.push_back(Loj);
		m_ChainHSPLens.push_back(Len);

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
#if TRACE_XDROP
	{
	uint M = SIZE(m_MuKmerHSPLois);
	Log("m_ChainCount=%u\n", M);
	for (uint Idx = 0; Idx < M; ++Idx)
		LogHSP(Idx);
	}
#endif
	}

void MuKmerFilter::Stats()
	{
	Log("MKF %u pairs, %u with HSPs (%.1f%%), box alns %u fract %.3g",
				m_PairCount, m_PairWithHSPCount,
				GetPct(m_PairWithHSPCount,m_PairCount),
				m_BoxAlnCount, m_SumBoxFract/(m_BoxAlnCount+1));
	Log(", hsps/chain %.1f",
				double(m_ChainHSPCount)/double(m_ChainCount+0.1));
	Log("\n");
	}

void MuKmerFilter::SetBagQ(const ChainBag &BagQ)
	{
	if (m_ptrBagQ == &BagQ)
		{
		assert(m_LabelQ == BagQ.m_ptrChain->m_Label);
		return;
		}

	m_ptrBagQ = &BagQ;
	m_LabelQ = BagQ.m_ptrChain->m_Label;

	m_ptrMuLettersQ = BagQ.m_ptrMuLetters;
	m_ptrMuKmersQ = BagQ.m_ptrMuKmers;
	if (m_ptrKmerHashTableQ == 0)
		m_ptrKmerHashTableQ = myalloc(uint16_t, m_DictSize*HASHW);
	memcpy(m_ptrKmerHashTableQ, 
		   BagQ.m_ptrKmerHashTableQ,
		   m_DictSize*HASHW*sizeof(uint16_t));

	asserta(m_ptrMuLettersQ != 0);
	asserta(m_ptrMuKmersQ != 0);
	asserta(m_ptrKmerHashTableQ != 0);
	}

void MuKmerFilter::LogHSP(uint Idx) const
	{
	asserta(Idx < SIZE(m_MuKmerHSPLois));
	int Loi = m_MuKmerHSPLois[Idx];
	int Loj = m_MuKmerHSPLojs[Idx];
	int Len = m_MuKmerHSPLens[Idx];
	int Score = m_MuKmerHSPScores[Idx];
	const byte *MuLettersQ = m_ptrMuLettersQ->data();
	const byte *MuLettersT = m_ptrMuLettersT->data();

	const char *SeqQ = m_DA->m_ChainA->m_Seq.c_str();
	const char *SeqT = m_DA->m_ChainB->m_Seq.c_str();

	int Score2 = 0;
	for (int Col = 0; Col < Len; ++Col)
		{
		byte q = MuLettersQ[Loi+Col];
		byte t = MuLettersT[Loj+Col];
		Score2 += IntScoreMx_Mu[q][t];
		}
	asserta(Score2 == Score);

	Log("\nMuKmerFilter::LogHSP(%u) Loi=%u Loj=%u Len=%u Score=%d\n",
		Idx, Loi, Loj, Len, Score);

	Log("Q/mu ");
	for (int Col = 0; Col < Len; ++Col)
		Log("%c", g_LetterToCharMu[MuLettersQ[Loi+Col]]);
	Log("\n");
	Log("T/mu ");
	for (int Col = 0; Col < Len; ++Col)
		Log("%c", g_LetterToCharMu[MuLettersT[Loj+Col]]);
	Log("\n");

	Log("Q/aa ");
	for (int Col = 0; Col < Len; ++Col)
		Log("%c", SeqQ[Loi+Col]);
	Log("\n");
	Log("T/aa ");
	for (int Col = 0; Col < Len; ++Col)
		Log("%c", SeqT[Loj+Col]);
	Log("\n");
	}
