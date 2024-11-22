#include "myutils.h"
#include "diagtrainer.h"
#include "diag.h"
#include "quarts.h"
#include "binner.h"

static void OnPairT(Trainer &Tr, uint PairIdx)
	{
	DiagTrainer &DT = (DiagTrainer &) Tr;
	uint ChainIdxQ = DT.m_ChainIdxsQ[PairIdx];
	uint ChainIdxR = DT.m_ChainIdxsR[PairIdx];
	DT.OnPair(ChainIdxQ, ChainIdxR, true);
	}

const byte *DiagTrainer::GetMuSeq(uint ChainIdx, uint &L)
	{
	const PDBChain &Chain = *m_Chains[ChainIdx];
	m_D.Init(Chain);
	vector<uint> Letters;
	m_D.GetMuLetters(Letters);
	L = Chain.GetSeqLength();
	byte *s = myalloc(byte, L);
	for (uint i = 0; i < L; ++i)
		{
		uint Letter = Letters[i];
		asserta(Letter < 36);
		s[i] = byte(Letter);
		}
	return s;
	}

static void OnPairF(Trainer &Tr, uint ChainIdxQ, uint ChainIdxR)
	{
	DiagTrainer &DT = (DiagTrainer &) Tr;
	DT.OnPair(ChainIdxQ, ChainIdxR, false);
	}

void DiagTrainer::Init(const DSSParams &Params, const string &PairAlnFN,
	const string &ChainsFN)
	{
	m_Params = &Params;
	m_D.SetParams(*m_Params);

	Trainer::Init(PairAlnFN, ChainsFN);
	}

bool DiagTrainer::SeedMatch(const byte *Q, int i, uint LQ,
							const byte *R, int j, uint LR) const
	{
	return false;
	}

bool DiagTrainer::GetKmer(const byte *Q, int i, uint QL, uint *Letters) const
	{
	if (i + m_k > QL)
		return false;
	for (uint offset = 0; offset < m_k; ++offset)
		Letters[offset] = Q[i+offset];
	return true;
	}

void DiagTrainer::OnPair(uint ChainIdxQ, uint ChainIdxR, bool IsT)
	{
	IsT ? ++m_NT : ++m_NF;
	
	uint LQ, LR;
	const byte *Q = GetMuSeq(ChainIdxQ, LQ);
	const byte *R = GetMuSeq(ChainIdxR, LR);

	m_DH.SetQ(Q, LQ);
	m_DH.SetT(R, LR);

	diag dg(LQ, LR);
	int mind = dg.getmind();
	int maxd = dg.getmaxd();
	int MaxScore = 0;
	int MaxDiag = 0;
	for (int d = mind; d <= maxd; ++d)
		{
		int Lo, Len;
		int Score = m_DH.Search(d, Lo, Len);
		if (Score > MaxScore)
			{
			MaxScore = Score;
			MaxDiag = d;
			}
		}

	int i = dg.getmini(MaxDiag);
	int j = dg.getminj(MaxDiag);
	int len = dg.getlen(MaxDiag);

	uint KmerLettersQ[16];
	uint KmerLettersR[16];
	uint SeedHitCount = 0;
	for (int offset = 0; offset < len; ++offset)
		{
		if (!GetKmer(Q, i+offset, LQ, KmerLettersQ))
			continue;
		if (!GetKmer(R, j+offset, LR, KmerLettersR))
			continue;

		extern const short * const *Mu_S_ij_short;
		int Score = 0;
		for (uint KmerPos = 0; KmerPos < m_k; ++KmerPos)
			{
			uint q = KmerLettersQ[KmerPos];
			uint r = KmerLettersR[KmerPos];
			assert(q < 36);
			assert(r < 36);
			Score += Mu_S_ij_short[q][r];
			if (Score > m_MinKmerScore)
				++SeedHitCount;
			}
		}

	if (SeedHitCount >= m_MinSeedHits && MaxScore >= m_MinDiagScore)
		{
		if (IsT)
			++m_SeedCountT;
		else
			++m_SeedCountF;
		}

	if (IsT)
		m_PosScores.push_back((float) MaxScore);
	else
		m_NegScores.push_back((float) MaxScore);
	}

void cmd_train_diags_scoredist()
	{
	asserta(optset_output);
	opt_fast = true;
	optset_fast = true;

	DSSParams Params;
	Params.SetFromCmdLine(10000);

	DiagTrainer DT;
	DT.Init(Params, g_Arg1, opt_train_cal);
	DT.EnumChainPairsT(OnPairT);
	DT.EnumChainPairsF(Scop40_IsTP_SF, OnPairF);

	ProgressLog("Seeds T %7u (%.1f%%)\n", DT.m_SeedCountT, GetPct(DT.m_SeedCountT, DT.m_NT));
	ProgressLog("Seeds F %7u (%.1f%%)\n", DT.m_SeedCountF, GetPct(DT.m_SeedCountF, DT.m_NF));

	QuartsFloat QuT, QuF;
	GetQuartsFloat(DT.m_PosScores, QuT);
	GetQuartsFloat(DT.m_NegScores, QuF);
	QuT.ProgressLogMe();
	QuF.ProgressLogMe();

	ProgressLog("True bins...");
	Binner <float>BT(DT.m_PosScores, 30, 0, 300);
	ProgressLog("\nFalse bins...");
	Binner <float>BF(DT.m_NegScores, 30, 0, 300);

	ProgressLog("\nWriting %s\n", opt_output);
	FILE *f = CreateStdioFile(opt_output);
	BT.ToTsv(f, "T");
	BF.ToTsv(f, "F");
	CloseStdioFile(f);
	}
