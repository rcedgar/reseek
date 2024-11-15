#include "myutils.h"
#include "seedtrainer.h"
#include "scop40bench.h"
#include "quarts.h"
#include "binner.h"

static inline bool isgap(char c) { return c == '-' || c == '.'; }

static inline bool OnSameDiagonal(uint MaxDist, uint PosQ1, uint PosR1,
  uint PosQ2, uint PosR2)
	{
	int Diag1 = int(PosQ1) - int(PosR1);
	int Diag2 = int(PosQ2) - int(PosR2);
	if (Diag1 != Diag2)
		return false;
	int Dist = abs(int(PosQ1) - int(PosQ2));
	if (Dist > int(MaxDist))
		return false;
	return true;
	}

static void OnPairT(Trainer &Tr, uint PairIdx)
	{
	SeedTrainer &ST = (SeedTrainer &) Tr;
	ST.OnPairT(PairIdx);
	}

static void OnPairF(Trainer &Tr, uint ChainIdxQ, uint ChainIdxR)
	{
	SeedTrainer &ST = (SeedTrainer &) Tr;
	ST.OnPairF(ChainIdxQ, ChainIdxR);
	}

float SeedTrainer::GetAaScore(byte Letter1, byte Letter2) const
	{
	return g_ScoreMxs2[FEATURE_AA][Letter1][Letter2];
	}

float SeedTrainer::GetMuScore(byte Letter1, byte Letter2) const
	{
	extern float ScoreMx_Mu[36][36];
	return ScoreMx_Mu[Letter1][Letter2];
	}

void SeedTrainer::GetAlignedKmerStarts(uint PairIdx,
	vector<uint> &PosQs, vector<uint> &PosRs) const
	{
	PosQs.clear();
	PosRs.clear();

	const string &RowQ = m_RowsQ[PairIdx];
	const string &RowR = m_RowsR[PairIdx];
	uint ColCount = SIZE(RowQ);

	uint ChainIdxQ = m_ChainIdxsQ[PairIdx];
	uint ChainIdxR = m_ChainIdxsR[PairIdx];

	uint LQ = GetSeqLength(ChainIdxQ);
	uint LR = GetSeqLength(ChainIdxR);

#if DEBUG
	{
	uint nq = 0;
	uint nr = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		if (isalpha(RowQ[Col])) ++nq;
		if (isalpha(RowR[Col])) ++nr;
		}
	assert(nq == LQ);
	assert(nr == LR);
	}
#endif

	string sq;
	string sr;
	asserta(SIZE(RowR) == ColCount);
	const string &Pattern = m_Params->m_PatternStr;
	uint PL = SIZE(Pattern);
	if (ColCount <= PL)
		return;
	uint PosQ = 0;
	uint PosR = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		bool Ok = true;
		for (uint i = 0; i < PL; ++i)
			{
			if (Col + i >= ColCount)
				{
				Ok = false;
				break;
				}
			char q = RowQ[Col+i];
			char r = RowR[Col+i];
			if (isgap(q) || isgap(r) || islower(q) || islower(r))
				{
				Ok = false;
				break;
				}
			}
		if (Ok)
			{
			assert(PosQ +  PL <= LQ);
			assert(PosR +  PL <= LR);
			PosQs.push_back(PosQ);
			PosRs.push_back(PosR);
			}
		if (!isgap(RowQ[Col])) ++PosQ;
		if (!isgap(RowR[Col])) ++PosR;
		}

	asserta(PosQ == LQ);
	asserta(PosR == LR);
	}

float SeedTrainer::GetAaKmerScore(uint Kmer1, uint Kmer2) const
	{
	float Score = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		byte Letter1 = Kmer1%20;
		byte Letter2 = Kmer2%20;
		Score += GetAaScore(Letter1, Letter2);
		Kmer1 = (Kmer1 - Letter1)/20;
		Kmer2 = (Kmer2 - Letter2)/20;
		}
	return Score;
	}

float SeedTrainer::GetMuKmerScore(uint Kmer1, uint Kmer2) const
	{
	float Score = 0;
	for (uint i = 0; i < m_k; ++i)
		{
		byte Letter1 = Kmer1%36;
		byte Letter2 = Kmer2%36;
		Score += GetMuScore(Letter1, Letter2);
		Kmer1 = (Kmer1 - Letter1)/36;
		Kmer2 = (Kmer2 - Letter2)/36;
		}
	return Score;
	}

bool SeedTrainer::IsPairSameDiag(
	const vector<uint> &PosQs1, const vector<uint> &PosRs1,
	const vector<uint> &PosQs2, const vector<uint> &PosRs2) const
	{
	for (uint qi1 = 0; qi1 < SIZE(PosQs1); ++qi1)
		{
		uint PosQi1 = PosQs1[qi1];
		for (uint ri1 = 0; ri1 < SIZE(PosRs1); ++ri1)
			{
			uint PosRi1 = PosRs1[ri1];
			int Diag1 = int(PosQi1) - int(PosRi1);
			for (uint qi2 = qi1+1; qi2 < SIZE(PosQs2); ++qi2)
				{
				uint PosQi2 = PosQs2[qi2];
				for (uint ri2 = ri1+1; ri2 < SIZE(PosRs2); ++ri2)
					{
					uint PosRi2 = PosRs2[ri2];
					int Diag2 = int(PosQi2) - int(PosRi2);
					if (Diag2 == Diag1)
						{
						int Dist = abs(int(PosQi2) - int(PosQi1));
						if (Dist >= int(m_k) && Dist <= int(m_MaxDiagDist))
							return true;
						}
					}
				}
			}
		}
	return false;
	}

/***
"score needed" = minimum of the two scores in a hit pair
MinScore
	= Lowest score threshold which finds one pair
	For TPs want threshold to be <= MinScore

MaxScore
	= Highest lower score over all pairs
	For FPs want threshold to be > MaxScore
***/
void SeedTrainer::FindIdenticalMuKmerPairsSameDiag(
  uint ChainIdxQ, uint ChainIdxR, float &MinScore, float &MaxScore) const
	{
	MinScore = FLT_MAX;
	MaxScore = FLT_MAX;
	const vector<uint> &MuKmersQ = m_MuKmersVec[ChainIdxQ];
	const vector<uint> &MuKmersR = m_MuKmersVec[ChainIdxR];

	map<uint, vector<uint> > KmerToPosVecQ;
	for (uint PosQ = 0; PosQ < SIZE(MuKmersQ); ++PosQ)
		{
		uint Kmer = MuKmersQ[PosQ];
		if (KmerToPosVecQ.find(Kmer) == KmerToPosVecQ.end())
			{
			vector<uint> v;
			KmerToPosVecQ[Kmer] = v;
			}
		KmerToPosVecQ[Kmer].push_back(PosQ);
		}

	vector<uint> IdenticalKmers;
	map<uint, vector<uint> > KmerToPosVecR;
	for (uint PosR = 0; PosR < SIZE(MuKmersR); ++PosR)
		{
		uint Kmer = MuKmersR[PosR];
		if (KmerToPosVecQ.find(Kmer) != KmerToPosVecQ.end())
			{
			IdenticalKmers.push_back(Kmer);
			if (KmerToPosVecR.find(Kmer) == KmerToPosVecR.end())
				{
				vector<uint> v;
				KmerToPosVecR[Kmer] = v;
				}
			KmerToPosVecR[Kmer].push_back(PosR);
			}
		}

	const uint N = SIZE(IdenticalKmers);
	uint PairSameDiagCount = 0;
	uint Count = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint Kmer_i = IdenticalKmers[i];
		const vector<uint> PosVecQ_i = KmerToPosVecQ[Kmer_i];
		const vector<uint> PosVecR_i = KmerToPosVecR[Kmer_i];
		for (uint j = i+1; j < N; ++j)
			{
			uint Kmer_j = IdenticalKmers[j];
			const vector<uint> PosVecQ_j = KmerToPosVecQ[Kmer_j];
			const vector<uint> PosVecR_j = KmerToPosVecR[Kmer_j];
			bool Yes = IsPairSameDiag(PosVecQ_i, PosVecR_i,
				PosVecQ_j, PosVecR_j);
			if (Yes)
				{
				float Score_i = GetMuKmerScore(Kmer_i, Kmer_i);
				float Score_j = GetMuKmerScore(Kmer_j, Kmer_j);
				float MinScore_ij = min(Score_i, Score_j);
				if (MinScore == FLT_MAX)
					{
					asserta(MaxScore == FLT_MAX);
					MinScore = MinScore_ij;
					MaxScore = MinScore_ij;
					}
				else
					{
					MinScore = min(MinScore, MinScore_ij);
					MaxScore = max(MaxScore, MinScore_ij);
					}
				}
			}
		}
	}

void SeedTrainer::OnPairT(uint PairIdx)
	{
	++m_NT;

	uint ChainIdxQ = m_ChainIdxsQ[PairIdx];
	uint ChainIdxR = m_ChainIdxsR[PairIdx];

	const float w = m_AaWeight;

	const vector<uint> &AaKmersQ = m_AaKmersVec[ChainIdxQ];
	const vector<uint> &AaKmersR = m_AaKmersVec[ChainIdxR];

	const vector<uint> &MuKmersQ = m_MuKmersVec[ChainIdxQ];
	const vector<uint> &MuKmersR = m_MuKmersVec[ChainIdxR];

	vector<uint> PosQs, PosRs;
	GetAlignedKmerStarts(PairIdx, PosQs, PosRs);
	const uint n = SIZE(PosQs);
	asserta(SIZE(PosRs) == n);

	uint KmerCount = min(SIZE(AaKmersQ), SIZE(AaKmersR));
	for (uint i = 0; i < n; ++i)
		{
		uint PosQ = PosQs[i];
		uint PosR = PosRs[i];
		float AaKmerScore = GetAaKmerScore(AaKmersQ[PosQ], AaKmersR[PosR]);
		float MuKmerScore = GetMuKmerScore(MuKmersQ[PosQ], MuKmersR[PosR]);
		float Score = w*AaKmerScore + MuKmerScore;
		m_PosScores.push_back(Score);
		}

	float MinScore, MaxScore_NOTUSED;
	 FindIdenticalMuKmerPairsSameDiag(
		 ChainIdxQ, ChainIdxR, MinScore, MaxScore_NOTUSED);
	if (MinScore != FLT_MAX)
		m_TwoHitIdenticalMuKmer_MinScores_T.push_back(MinScore);
	}

void SeedTrainer::OnPairF(uint ChainIdxQ, uint ChainIdxR)
	{
	++m_NF;
	const float w = m_AaWeight;

	const vector<uint> &AaKmersQ = m_AaKmersVec[ChainIdxQ];
	const vector<uint> &AaKmersR = m_AaKmersVec[ChainIdxR];

	const vector<uint> &MuKmersQ = m_MuKmersVec[ChainIdxQ];
	const vector<uint> &MuKmersR = m_MuKmersVec[ChainIdxR];

	uint KmerCount = min(SIZE(AaKmersQ), SIZE(AaKmersR));
	for (uint i = 0; i < KmerCount; ++i)
		{
		uint j = randu32()%KmerCount;
		float AaKmerScore = GetAaKmerScore(AaKmersQ[i], AaKmersR[j]);
		float MuKmerScore = GetMuKmerScore(MuKmersQ[i], MuKmersR[j]);
		float Score = w*AaKmerScore + MuKmerScore;
		m_NegScores.push_back(Score);
		}

	float MinScore_NOTUSED, MaxScore;
	 FindIdenticalMuKmerPairsSameDiag(
		 ChainIdxQ, ChainIdxR, MinScore_NOTUSED, MaxScore);
	if (MaxScore != FLT_MAX)
		m_TwoHitIdenticalMuKmer_MaxScores_F.push_back(MaxScore);
	}

void cmd_train_seeds_scoredist()
	{
	asserta(optset_output);
	opt_fast = true;
	optset_fast = true;

	DSSParams Params;
	Params.SetFromCmdLine(10000);

	SeedTrainer ST;
	ST.m_AaWeight = 1;
	if (optset_aaweight) ST.m_AaWeight = (float) opt_aaweight;
	ST.Init(Params, g_Arg1, opt_train_cal);
	ST.EnumChainPairsT(OnPairT);
	ST.EnumChainPairsF(Scop40_IsTP_SF, OnPairF);

	QuartsFloat QT;
	QuartsFloat QF;
	GetQuartsFloat(ST.m_PosScores, QT);
	GetQuartsFloat(ST.m_NegScores, QF);
	QT.LogMe();
	QF.LogMe();
	Log("NT %u, NF %u\n", ST.m_NT, ST.m_NF);

	const uint BINS = 33;
	Binner<float> BT(ST.m_PosScores, BINS, -5, +5);
	Binner<float> BF(ST.m_NegScores, BINS, -5, +5);

	uint NT = BT.GetTotalCount();
	uint NF = BF.GetTotalCount();
	const vector<uint> &bts = BT.GetBins();
	const vector<uint> &bfs = BF.GetBins();

	FILE *f = CreateStdioFile(opt_output);
	fprintf(f, "Kmer\tBinLo\tNT\tNF\tFreqT\tFreqF\n");
	for (uint b = 1; b + 1 < BINS; ++b)
		{
		float Lo = BT.GetBinLo(b);
		uint nt = bts[b];
		uint nf = bfs[b];
		float freqt = float(nt)/NT;
		float freqf = float(nf)/NF;
		fprintf(f, "Kmer\t%.3f\t%u\t%u\t%.3g\t%.3g\n",
			Lo, nt, nf, freqt, freqf);
		}

	Binner<float> BT2(ST.m_TwoHitIdenticalMuKmer_MinScores_T, BINS, -5, +15);
	Binner<float> BF2(ST.m_TwoHitIdenticalMuKmer_MaxScores_F, BINS, -5, +15);
	vector<uint> bts2;
	vector<uint> bfs2;
	BT2.GetAccumBinsReverse(bts2);
	BF2.GetAccumBinsReverse(bfs2);
	fprintf(f, "IdMu2\tBinLo\tNT\tNF\tFreqT\tFreqF\n");
	for (uint b = 1; b + 1 < BINS; ++b)
		{
		float Lo = BT2.GetBinLo(b);
		uint nt = bts2[b];
		uint nf = bfs2[b];
		float freqt = float(nt)/NT;
		float freqf = float(nf)/NF;
		fprintf(f, "IdMu2\t%.3f\t%u\t%u\t%.3g\t%.3g\n",
			Lo, nt, nf, freqt, freqf);
		}

	CloseStdioFile(f);
	}
