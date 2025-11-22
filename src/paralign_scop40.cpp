#include "myutils.h"
#include "paralign.h"
#include "seqdb.h"
#include "alpha.h"
#include "scop40bench.h"
#include "triangle.h"

#define POKE_MU	1

void cmd_paralign_scop40()
	{
	asserta(optset_bins);
	asserta(optset_lookup);
	SCOP40Bench SB;
	SB.ReadLookup(opt(lookup));

	const string &FastaFN = g_Arg1;
	SeqDB Seqs;
	Seqs.FromFasta(FastaFN);

	if (opt(bins) == 8)
		Paralign::m_Bits = 8;
	else if (opt(bins) == 16)
		Paralign::m_Bits = 16;
	else
		Die("-bins must be 8 or 16");

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

#if POKE_MU
	{
	extern int8_t IntScoreMx_Mu[36][36];
	vector<vector<int> > ScoreMx(36);
	for (uint i = 0; i < 36; ++i)
		{
		ScoreMx[i].resize(36);
		for (uint j = 0; j < 36; ++j)
			ScoreMx[i][j] = IntScoreMx_Mu[i][j];
		}
	Paralign::SetMatrix(ScoreMx, 2, 1, 777);
	}
#else
	Paralign::SetMu();
#endif
	void Log_parasail_mu_matrix(const parasail_matrix_t &mx);//@@TODO
	Log_parasail_mu_matrix(Paralign::m_matrix);

	uint PairCount = SeqCount*(SeqCount-1)/2 + SeqCount;
	uint PairCount2 = triangle_get_k(SeqCount) + 1;
	ProgressLog("PairCount %u %u\n", PairCount, PairCount2);
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
	}
