#include "myutils.h"
#include "paralign.h"
#include "seqdb.h"
#include "alpha.h"
#include "scop40bench.h"
#include "triangle.h"

void FixMuByteSeq(vector<byte> &ByteSeq);

static void GetByteSeqs_Mux(const string &FN,
	vector<string> &Labels,
	vector<vector<byte> > &ByteSeqs)
	{
	//vector<FEATURE> Fs;
	//Fs.push_back(FEATURE_SS3);
	//Fs.push_back(FEATURE_NENSS3);
	//Fs.push_back(FEATURE_RENDist4);

	vector<PDBChain *> Chains;
	ReadChains(FN, Chains);
	const uint ChainCount = SIZE(Chains);
	DSS D;
	ByteSeqs.clear();
	ByteSeqs.resize(ChainCount);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		Labels.push_back(Chain.m_Label);
		D.Init(Chain);
		vector<byte> &ByteSeq = ByteSeqs[ChainIdx];
		ByteSeq.reserve(L);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint Letter_SS3 = D.GetFeature(FEATURE_SS3, Pos);
			uint Letter_NENSS3 = D.GetFeature(FEATURE_NENSS3, Pos);
			uint Letter_RENDist4 = D.GetFeature(FEATURE_RENDist4, Pos);
			byte Letter = Letter_SS3 + Letter_NENSS3*3 + Letter_RENDist4*3*3;
			byte Letter2 = D.Get_Mu(Pos);
			assert(Letter == Letter2);
			asserta(Letter < 36);
			ByteSeq.push_back(Letter);
			}
		if (opt(fixmubyteseq))
			FixMuByteSeq(ByteSeq);
		}
	}

static void GetByteSeqs_Mu(const string &FN,
	vector<string> &Labels,
	vector<vector<byte> > &ByteSeqs)
	{
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
		if (opt(fixmubyteseq))
			FixMuByteSeq(ByteSeq);
		}
	}

void cmd_paralign_scop40_floatmu()
	{
	asserta(!optset_alpha);
	asserta(optset_intopen);
	asserta(optset_intext);
	asserta(optset_lookup);
	SCOP40Bench SB;
	SB.ReadLookup(opt(lookup));

	const string &ChainsFN = g_Arg1;
	vector<string> Labels;
	vector<vector<byte> > ByteSeqs;
	GetByteSeqs_Mu(ChainsFN, Labels, ByteSeqs);

	const uint SeqCount = SIZE(ByteSeqs);

 	Paralign::SetMu_musubstmx();

	uint PairCount = SeqCount*(SeqCount-1)/2 + SeqCount;
	uint PairCount2 = triangle_get_K(SeqCount) + 1;
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
		PA.SetQueryNoProfile(Label_i, ByteSeq_i.data(), L_i);

		for (uint j = i+1; j < SeqCount; ++j)
			{
			const string &Label_j = Labels[j];
			const vector<byte> &ByteSeq_j = ByteSeqs[j];
			uint L_j = SIZE(ByteSeq_j);
			PA.Align_SWFast(Label_j, ByteSeq_j.data(), L_j);
#pragma omp critical
			{
			Label1s.push_back(Label_i);
			Label2s.push_back(Label_j);
			Scores.push_back(PA.m_SWFastScore);

			Label1s.push_back(Label_j);
			Label2s.push_back(Label_i);
			Scores.push_back(PA.m_SWFastScore);
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

	if (optset_output)
		{
		FILE *f = CreateStdioFile(opt(output));
		const uint HitCount = SB.GetHitCount();
		const vector<uint> &Order = SB.m_ScoreOrder;
		asserta(SIZE(Order) == HitCount);
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
			fprintf(f, "%.1f\t%s\t%s\n", Score, Label1.c_str(), Label2.c_str());
			}
		CloseStdioFile(f);
		}
	}
