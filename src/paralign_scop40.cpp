#include "myutils.h"
#include "paralign.h"
#include "seqdb.h"
#include "alpha.h"
#include "scop40bench.h"
#include "triangle.h"

static void GetByteSeqs(const string &FN,
	vector<string> &Labels,
	vector<vector<byte> > &ByteSeqs)
	{
	FILE *f = CreateStdioFile("mu.tmp");
	vector<PDBChain *> Chains;
	ReadChains(FN, Chains);
	const uint ChainCount = SIZE(Chains);
	DSS D;
	ByteSeqs.clear();
	ByteSeqs.resize(ChainCount);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *Chains[ChainIdx];
		Labels.push_back(Chain.m_Label);
		D.Init(Chain);
		vector<byte> &ByteSeq = ByteSeqs[ChainIdx];
		D.GetMuLetters(ByteSeq);

	/////////////////////////////////////////
	// Hack because of K<->L bug in alpha.cpp
	/////////////////////////////////////////
	// unsigned char g_CharToLetterMu[256] =
	// ...
	//	9  ,         // [ 74] 'J'
	//	11 ,         // [ 75] 'L'
	//	10 ,         // [ 76] 'K'
	//	12 ,         // [ 77] 'M'
	//unsigned char g_LetterToCharMu[256] =
	// ...
	//	'J',           // [9]
	//	'L',           // [10]
	//	'K',           // [11]
	//	'M',           // [12]
	/////////////////////////////////////////
		for (uint i = 0; i < SIZE(ByteSeq); ++i)
			{
			byte Letter = ByteSeq[i];
			if (Letter == 11)
				ByteSeq[i] = 10;
			else if (Letter == 10)
				ByteSeq[i] = 11;
			}
	/////////////////////////////////////////
		}
	CloseStdioFile(f);
	}

void cmd_paralign_scop40()
	{
	asserta(optset_lookup);
	SCOP40Bench SB;
	SB.ReadLookup(opt(lookup));

	const string &ChainsFN = g_Arg1;
	vector<string> Labels;
	vector<vector<byte> > ByteSeqs;
	GetByteSeqs(ChainsFN, Labels, ByteSeqs);

	const uint SeqCount = SIZE(ByteSeqs);

	Paralign::SetMu();

	uint PairCount = SeqCount*(SeqCount-1)/2 + SeqCount;
	uint PairCount2 = triangle_get_k(SeqCount) + 1;
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
		if (i < SeqCount)
			ProgressStep(i, SeqCount, "Aligning");
		}
		if (i >= SeqCount)
			break;

		const string &Label_i = Labels[i];
		const vector<byte> &ByteSeq_i = ByteSeqs[i];
		uint L_i = SIZE(ByteSeq_i);
		PA.SetQuery(Label_i, ByteSeq_i.data(), L_i);

		for (uint j = i+1; j < SeqCount; ++j)
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

#if 0
	{
	const uint HitCount = SB.GetHitCount();
	const vector<uint> &Order = SB.m_ScoreOrder;
	asserta(SIZE(Order) == HitCount);
	FILE *f = CreateStdioFile("paralign_scop40.hits");
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
