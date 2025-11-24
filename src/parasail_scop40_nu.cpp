#include "myutils.h"
#include "paralign.h"
#include "alpha.h"
#include "scop40bench.h"
#include "triangle.h"
#include "nu.h"

void FixMuByteSeq(vector<byte> &ByteSeq);

void cmd_paralign_scop40_nu()
	{
	asserta(optset_lookup);
	SCOP40Bench SB;
	SB.ReadLookup(opt(lookup));

	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	const uint ChainCount = SIZE(Chains);

	Nu A;
	A.SetMu();

	vector<vector<byte> > ByteSeqs(ChainCount);
	vector<string> Labels;
	for (uint i = 0; i < ChainCount; ++i)
		{
		ProgressStep(i, ChainCount, "Byte seqs");
		const PDBChain &Chain = *Chains[i];
		const uint L = Chain.GetSeqLength();
		vector<byte> &ByteSeq = ByteSeqs[i];
		A.GetLetters(Chain, ByteSeq);
		if (opt(fixmubyteseq))
			FixMuByteSeq(ByteSeq);
		Labels.push_back(Chain.m_Label);
		}

	vector<vector<float> > NuMxf;
	A.GetScoreMx(NuMxf);

	vector<vector<int> > NuMxi;
	A.FloatMxToIntMx(NuMxf, 3, NuMxi);
	A.MxToSrci(g_fLog, "int", "NuMxi", 3, NuMxi);

	Paralign::SetMatrix(NuMxi, 10, 3, 777);

	uint PairCount = ChainCount*(ChainCount-1)/2 + ChainCount;
	uint PairCount2 = triangle_get_k(ChainCount) + 1;
	asserta(PairCount == PairCount2);
	uint ScoreDiffs = 0;
	vector<string> Label1s;
	vector<string> Label2s;
	vector<float> Scores;
	Label1s.reserve(PairCount);
	Label2s.reserve(PairCount);
	Scores.reserve(PairCount);
	const uint ThreadCount = GetRequestedThreadCount();
	vector<Paralign *> PAs;
	for (uint i = 0; i < ThreadCount; ++i)
		PAs.push_back(new Paralign);

	uint Next_i = 0;
#pragma omp parallel num_threads(ThreadCount)
	while (true)
		{
		uint ThreadIndex = GetThreadIndex();
		Paralign &PA = *PAs[ThreadIndex];
		uint i;
#pragma omp critical
		{
		i = Next_i++;
		if (i < ChainCount)
			ProgressStep(i, ChainCount, "Aligning");
		}
		if (i >= ChainCount)
			break;

		const string &Label_i = Labels[i];
		const vector<byte> &ByteSeq_i = ByteSeqs[i];
		uint L_i = SIZE(ByteSeq_i);
		PA.SetQuery(Label_i, ByteSeq_i.data(), L_i);

		for (uint j = i+1; j < ChainCount; ++j)
			{
			const string &Label_j = Labels[j];
			const vector<byte> &ByteSeq_j = ByteSeqs[j];
			uint L_j = SIZE(ByteSeq_j);

			PA.Align_ScoreOnly(Label_j, ByteSeq_j.data(), L_j);

#pragma omp critical
			{
			Label1s.push_back(Label_i);
			Label2s.push_back(Label_j);
			Scores.push_back((float) PA.m_Score);
			}
			}
		}

	ProgressLog("%u long, %u saturated, %u 8-bit, %u 16-bit\n",
		Paralign::m_TooLongCount.load(),
		Paralign::m_SaturatedCount.load(),
		Paralign::m_Count8.load(),
		Paralign::m_Count16.load());

	SB.SetHits(Label1s, Label2s, Scores);
	SB.m_SBS = SBS_OtherAlgoScore;
	SB.SetScoreOrder();
	SB.WriteOutput();

#if 1
	{
	const uint HitCount = SB.GetHitCount();
	const vector<uint> &Order = SB.m_ScoreOrder;
	asserta(SIZE(Order) == HitCount);
	FILE *f = CreateStdioFile("paralign_scop40_nu.hits");
	for (uint k = 0; k < HitCount; ++k)
		{
		uint i = k;
		uint DomIdx1 = SB.m_DomIdx1s[i];
		uint DomIdx2 = SB.m_DomIdx2s[i];
		const string &Label1 = SB.m_Doms[DomIdx1];
		const string &Label2 = SB.m_Doms[DomIdx2];
		if (Label1 == Label2)
			continue;
		float Score = SB.m_Scores[i];
		fprintf(f, "%.0f\t%s\t%s\n", Score, Label1.c_str(), Label2.c_str());
		}
	}
#endif
	}
