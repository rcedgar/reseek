#include "myutils.h"
#include "paralign.h"
#include "seqdb.h"
#include "alpha.h"
#include "scop40bench.h"
#include "triangle.h"

/***
reseek v2.8.win64 [dac28ca]
C:\src\reseek\src\Release\reseek.exe -paralign_scop40_mu ../out/a.mu -lookup ../out/a.lookup -log paralign_scop40_mu.log -intopen 20 -intext 10 -bins 16 
Started Sat Nov 22 14:31:13 2025
0 long, 0 saturated, 0 8-bit, 8419356 16-bit

SEPQ0.1=0.105 SEPQ1=0.168 SEPQ10=0.263 Area0=0.429 Sum3=0.725 ../out/a.mu

Elapsed time 00:27
Max memory 1.4Gb
Finished Sat Nov 22 14:31:40 2025

With POKE_MU=1 (old 8-bit settings)
SEPQ0.1=0.066 SEPQ1=0.106 SEPQ10=0.215 Area0=0.283 Sum3=0.506
***/

void cmd_paralign_scop40_mu()
	{
	asserta(optset_lookup);
	SCOP40Bench SB;
	SB.ReadLookup(opt(lookup));

	const string &FastaFN = g_Arg1;
	SeqDB Seqs;
	Seqs.FromFasta(FastaFN);

	const uint SeqCount = Seqs.GetSeqCount();
	vector<vector<byte> > ByteSeqs(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		ProgressStep(i, SeqCount, "Byte seqs");
		const string &Seq = Seqs.GetSeq(i);
		const uint L = Seqs.GetSeqLength(i);
		vector<byte> &ByteSeq = ByteSeqs[i];
		ByteSeq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			char c = Seq[Pos];
			byte Letter = g_CharToLetterMu[c];
			ByteSeq.push_back(Letter);
			}
		}

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

		const string &Label_i = Seqs.GetLabel(i);
		const string &Seq_i = Seqs.GetSeq(i);
		const byte *ByteSeq_i = ByteSeqs[i].data();
		uint L_i = Seqs.GetSeqLength(i);
		PA.SetQuery(Label_i, ByteSeq_i, L_i);

		for (uint j = i+1; j < SeqCount; ++j)
			{
			const string &Label_j = Seqs.GetLabel(j);
			const string &Seq_j = Seqs.GetSeq(j);
			const byte *ByteSeq_j = ByteSeqs[j].data();
			uint L_j = Seqs.GetSeqLength(j);
			PA.Align_ScoreOnly(Label_j, ByteSeq_j, L_j);
#pragma omp critical
			{
			Label1s.push_back(Label_i);
			Label2s.push_back(Label_j);
			Scores.push_back((float) PA.m_Score);

			Label1s.push_back(Label_j);
			Label2s.push_back(Label_i);
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
	FILE *f = CreateStdioFile("paralign_scop40_mu.hits");
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
