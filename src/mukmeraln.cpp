#include "myutils.h"
#include "dbsearcher.h"
#include "timing.h"

#if MUKMERS

//static const uint DictSize = 36*36*36*36*36;
//static const string PatternStr = "11111";
//static const uint DictSize = 36*36*36*36;
//static const string PatternStr = "1111";
static const uint DictSize = 36*36*36;
static const string PatternStr = "111";
extern int8_t IntScoreMx_Mu[36][36];

void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);

static void GetColToPos(const string &Path, bool IsQ, uint Lo, uint L,
						vector<uint> &ColToPos)
	{
	const uint ColCount = SIZE(Path);
	uint Pos = Lo;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		bool IsGap = false;
		if (IsQ)
			{
			if (c == 'I')
				IsGap = true;
			}
		else
			{
			if (c == 'D')
				IsGap = true;
			}
		if (IsGap)
			ColToPos.push_back(UINT_MAX);
		else
			{
			asserta(Pos < L);
			ColToPos.push_back(Pos++);
			}
		}
	}

static void GetPosToPosVecs(const string &Path, uint LoQ, uint LQ, uint LoT, uint LT,
							vector<uint> &PosQToPosT,
							vector<uint> &PosTToPosQ)
	{
	PosQToPosT.clear();
	PosTToPosQ.clear();

	const uint ColCount = SIZE(Path);
	vector<uint> ColToPosQ;
	vector<uint> ColToPosT;
	GetColToPos(Path, true, LoQ, LQ, ColToPosQ);
	GetColToPos(Path, false, LoT, LT, ColToPosT);
	asserta(SIZE(ColToPosQ) == ColCount);
	asserta(SIZE(ColToPosT) == ColCount);

	PosQToPosT.resize(LQ, UINT_MAX);
	PosTToPosQ.resize(LT, UINT_MAX);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint PosQ = ColToPosQ[Col];
		uint PosT = ColToPosT[Col];
		if (PosQ != UINT_MAX)
			{
			asserta(PosQ < LQ);
			PosQToPosT[PosQ] = PosT;
			}
		if (PosT != UINT_MAX)
			{
			asserta(PosT < LT);
			PosTToPosQ[PosT] = PosQ;
			}
		}
	}

void DBSearcher::MuKmerCmpHSPPath(DSSAligner &DA)
	{
	const string &Path = DA.m_Path;
	if (Path.empty())
		return;
	const uint N = SIZE(m_MuKmerBestLens);
	asserta(N > 0);
	uint LoQ = DA.m_LoA;
	uint LoT = DA.m_LoB;
	uint LQ = DA.GetL(true);
	uint LT = DA.GetL(false);
	//uint M, D, I;
	//GetPathCounts(Path, M, D, I);
	//asserta(LoQ + M + D <= LQ);
	//asserta(LoT + M + I <= LT);
	//vector<uint> PosQToPosT;
	//vector<uint> PosTToPosQ;
	//GetPosToPosVecs(DA.m_Path, LoQ, LQ, LoT, LT, PosQToPosT, PosTToPosQ);

	//int i = m_MuKmerBestLoi;
	//int j = m_MuKmerBestLoj;
	//Log("\n");
	//Log(">%s, %s (%.3g)\n",
	//	DA.m_ChainA->m_Label.c_str(), DA.m_ChainB->m_Label.c_str(),
	//	DA.m_EvalueA);
	//for (int k = 0; k < m_MuKmerBestLen; ++k)
	//	{
	//	uint j2 = PosQToPosT[i];
	//	uint i2 = PosTToPosQ[j];
	//	Log("i=%5d  i2=%5d  j=%5d  j2=%5d\n", i, i2, j, j2);
	//	++i;
	//	++j;
	//	}
	Log("@CMP@");
	Log("\t%.3g", DA.m_EvalueA);
	Log("\t%u", N);
	Log("\t%u", LoQ);
	Log("\t%u", LQ);
	Log("\t%u", LoT);
	Log("\t%u", LT);
	Log("\t%s", Path.c_str());
	for (uint i = 0; i < N; ++i)
		{
		Log("\t%d", m_MuKmerBestLois[i]);
		Log("\t%d", m_MuKmerBestLojs[i]);
		Log("\t%d", m_MuKmerBestLens[i]);
		}
	Log("\n");
	}

int DBSearcher::MuXDrop(int PosQ, int LQ, int PosT, int LT, int X,
						int &Loi, int &Loj, int &Len)
	{
	StartTimer(MuXDrop);
	Loi = PosQ;
	Loj = PosT;
	Len = 0;
	int i = PosQ;
	int j = PosT;
	const byte *Q = m_MuLettersQ.data();
	const byte *T = (*m_ptrMuLettersT).data();

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
	++m_MuKmerFilterPairCount;
	const uint KmerCountT = SIZE(MuKmersT);
	int LQ = int(m_ChainQ->GetSeqLength());
	int LT = int(m_ChainT->GetSeqLength());
	bool FoundHSP = false;
	int BestHSPScore = 0;
	m_MuKmerBestLois.clear();
	m_MuKmerBestLojs.clear();
	m_MuKmerBestLens.clear();
	for (uint PosT = 0; PosT < KmerCountT; ++PosT)
		{
		uint KmerT = MuKmersT[PosT];
		assert(KmerT < DictSize);
		uint PosQ = m_KmerHashTableQ[KmerT];
		if (PosQ != 0xffff)
			{
			EndTimer(MuKmerAln);
			int Loi, Loj, Len;
			int Score = MuXDrop(int(PosQ), LQ, int(PosT), LT, 10, Loi, Loj, Len);
			StartTimer(MuKmerAln);
			if (Score >= 75)
				{
				if (Score > BestHSPScore)
					{
					FoundHSP = true;
					bool OldHSP = false;
					for (uint i = 0; i < SIZE(m_MuKmerBestLois); ++i)
						{
						if (Loi == m_MuKmerBestLois[i])
							{
							OldHSP = true;
							break;
							}
						}
					if (!OldHSP)
						{
						m_MuKmerBestLois.push_back(Loi);
						m_MuKmerBestLojs.push_back(Loj);
						m_MuKmerBestLens.push_back(Len);
						}
					}
				}
			}
		}
	EndTimer(MuKmerAln);
	if (FoundHSP)
		++m_MuKmerFilterHitCount;
	return FoundHSP;
	}
#endif // MUKMERS