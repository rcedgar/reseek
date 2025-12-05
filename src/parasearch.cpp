#include "myutils.h"
#include "parasearch.h"
#include "triangle.h"
#include "nu.h"

/////////////////////////////////////////
// Hack because of K<->L bug in alpha.cpp
//   baked into mumx_data.cpp
/////////////////////////////////////////
// unsigned char g_CharToLetterMu[256] =
// ...
//	9  ,         // [ 74] 'J'
//	11 ,         // [ 75] 'L'
//	10 ,         // [ 76] 'K'
//	12 ,         // [ 77] 'M'
//
//unsigned char g_LetterToCharMu[256] =
// ...
//	'J',           // [9]
//	'L',           // [10]
//	'K',           // [11]
//	'M',           // [12]
/////////////////////////////////////////
void FixMuByteSeq(vector<byte> &ByteSeq)
	{
	for (uint i = 0; i < SIZE(ByteSeq); ++i)
		{
		byte Letter = ByteSeq[i];
		if (Letter == 11)
			ByteSeq[i] = 10;
		else if (Letter == 10)
			ByteSeq[i] = 11;
		}
	}

vector<FEATURE> ParaSearch::m_NuFs;

void ParaSearch::AppendHit(uint i, uint j, float Score)
	{
	uint k = triangle_ij_to_k(i, j, m_SeqCount);
	m_Scores[k] = Score;
	}

void ParaSearch::Align(uint ThreadIdx, uint i, uint j)
	{
	const string &Label_j = m_Labels[j];
	const vector<byte> &ByteSeq_j = m_ByteSeqs[j];
	uint L_j = SIZE(ByteSeq_j);

	if (i != m_QueryIdxs[ThreadIdx])
		{
		SetQuery(ThreadIdx, i);
		m_QueryIdxs[ThreadIdx] = i;
		}

	Paralign &PA = *m_PAs[ThreadIdx];
	PA.m_LabelT = Label_j;
	PA.m_T = ByteSeq_j.data();
	PA.m_LT = L_j;

	if (m_AlignMethod == "sw")
		{
		PA.Align_SWFast(Label_j, ByteSeq_j.data(), L_j);
		AppendHit(i, j, PA.m_SWFastScore);
		}
	else if (m_AlignMethod == "para")
		{
		PA.Align_ScoreOnly(Label_j, ByteSeq_j.data(), L_j);
		AppendHit(i, j, (float) PA.m_Score);
		}
	else
		Die("m_AlignMethod=%s", m_AlignMethod.c_str());
	}

void ParaSearch::SetQuery(uint ThreadIdx, uint i)
	{
	const string &Label_i = m_Labels[i];
	const vector<byte> &ByteSeq_i = m_ByteSeqs[i];
	uint L_i = SIZE(ByteSeq_i);
	Paralign &PA = *m_PAs[ThreadIdx];
	PA.m_Q = ByteSeq_i.data();
	PA.m_LQ = L_i;

	if (m_AlignMethod == "para")
		PA.SetQueryProfile(Label_i, ByteSeq_i.data(), L_i);
	else if (m_AlignMethod == "sw")
		;
	else
		Die("m_AlignMethod=%s", m_AlignMethod.c_str());
	}

void ParaSearch::Search(const string &AlignMethod)
	{
	m_AlignMethod = AlignMethod;

	ProgressLog("Search %s %s %s\n",
		m_AlignMethod.c_str(),
		m_SubstMxName.c_str(),
		m_ByteSeqMethod.c_str());

	uint PairCount2 = triangle_get_K(m_SeqCount) + 1;
	asserta(m_PairCount == PairCount2);
	m_Scores = myalloc(float, m_PairCount);
	const uint ThreadCount = GetRequestedThreadCount();
	m_PAs.clear();
	m_QueryIdxs.clear();
	for (uint i = 0; i < ThreadCount; ++i)
		{
		Paralign *PA = new Paralign;
		m_PAs.push_back(PA);
		m_QueryIdxs.push_back(UINT_MAX);
		}

	atomic<uint> Counter = 0;
	ProgressStep(0, m_PairCount, "Aligning");

#pragma omp parallel num_threads(ThreadCount)
	{
	uint ThreadIdx = GetThreadIndex();
#pragma omp for
	for (int PairIdx = 0; PairIdx < int(m_PairCount); ++PairIdx)
		{
		++Counter;
		if (Counter.load()%1000 == 0)
			{
#pragma omp critical
				{
				ProgressStep(Counter.load(), m_PairCount, "Aligning");
				}
			}

		uint i, j;
		triangle_k_to_ij(PairIdx, m_SeqCount, i, j);
		Align(ThreadIdx, i, j);
		}
	}

	ProgressStep(m_PairCount-1, m_PairCount, "Aligning");
	ProgressLog("%u saturated, %u 8-bit, %u 16-bit, %u SW\n",
		Paralign::m_SaturatedCount.load(),
		Paralign::m_Count8.load(),
		Paralign::m_Count16.load(),
		Paralign::m_CountSWFast.load());
	}

void ParaSearch::GetByteSeqs(const string &FN, const string &Method)
	{
	m_ByteSeqMethod = Method;
	if (Method == "muletters")
		GetByteSeqs_muletters(FN);
	else if (Method == "dss3")
		GetByteSeqs_dss3(FN);
	else if (Method == "numu")
		GetByteSeqs_numu(FN);
	else if (Method == "nuletters")
		GetByteSeqs_nu(FN);
	else
		Die("GetByteSeqs(%s)", Method.c_str());

	m_SeqCount = SIZE(m_ByteSeqs);
	m_PairCount = m_SeqCount*(m_SeqCount-1)/2 + m_SeqCount;
	}

// Construct Mu from components to validate that it
// reproduces DSS::GetMu(). Otherwise this is redundant,
// better to use ParaSearch::GetByteSeqs_DSS().
void ParaSearch::GetByteSeqs_dss3(const string &FN)
	{
	ReadChains(FN, m_Chains);
	const uint ChainCount = SIZE(m_Chains);

	DSS D;
	m_ByteSeqs.clear();
	m_ByteSeqs.resize(ChainCount);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *m_Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		m_Labels.push_back(Chain.m_Label);
		D.Init(Chain);
		vector<byte> &ByteSeq = m_ByteSeqs[ChainIdx];
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

// Use DSS::GetMuLetters()
void ParaSearch::GetByteSeqs_muletters(const string &FN)
	{
	ReadChains(FN, m_Chains);
	const uint ChainCount = SIZE(m_Chains);
	DSS D;
	m_ByteSeqs.clear();
	m_ByteSeqs.resize(ChainCount);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "DSS::GetMuLetters()");
		const PDBChain &Chain = *m_Chains[ChainIdx];
		m_Labels.push_back(Chain.m_Label);
		D.Init(Chain);
		vector<byte> &ByteSeq = m_ByteSeqs[ChainIdx];
		D.GetMuLetters(ByteSeq);
		if (opt(fixmubyteseq))
			FixMuByteSeq(ByteSeq);
		}
	}

void ParaSearch::GetByteSeqs_nu(const string &FN)
	{
	m_ByteSeqs.clear();
	m_Labels.clear();

	vector<float> Weights;
	const uint NF = SIZE(m_NuFs);
	for (uint i = 0; i < NF; ++i)
		Weights.push_back(float(1)/i);
	
	Nu TheNu;
	TheNu.SetComponents(m_NuFs, Weights);

	ReadChains(FN, m_Chains);
	const uint ChainCount = SIZE(m_Chains);
	m_ByteSeqs.resize(ChainCount);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "GetByteSeqs_nu()");
		const PDBChain &Chain = *m_Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		vector<byte> &ByteSeq = m_ByteSeqs[ChainIdx];
		TheNu.GetLetters(Chain, ByteSeq);
		if (opt(fixmubyteseq))
			FixMuByteSeq(ByteSeq);
		m_Labels.push_back(Chain.m_Label);
		}
	}

void ParaSearch::GetByteSeqs_numu(const string &FN)
	{
	m_ByteSeqs.clear();
	m_Labels.clear();

	Nu A;
	A.SetMu();

	ReadChains(FN, m_Chains);
	const uint ChainCount = SIZE(m_Chains);
	m_ByteSeqs.resize(ChainCount);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Nu::GetLetters(Mu)");
		const PDBChain &Chain = *m_Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		vector<byte> &ByteSeq = m_ByteSeqs[ChainIdx];
		A.GetLetters(Chain, ByteSeq);
		if (opt(fixmubyteseq))
			FixMuByteSeq(ByteSeq);
		m_Labels.push_back(Chain.m_Label);
		}
	}

void ParaSearch::Bench()
	{
	vector<string> Label1s;
	vector<string> Label2s;
	vector<float> Scores;
	Label1s.reserve(m_PairCount);
	Label2s.reserve(m_PairCount);
	Scores.reserve(m_PairCount);
	for (uint k = 0; k < m_PairCount; ++k)
		{
		uint i, j;
		triangle_k_to_ij(k, m_SeqCount, i, j);
		Label1s.push_back(m_Labels[i]);
		Label2s.push_back(m_Labels[j]);
		Scores.push_back(m_Scores[k]);

		Label1s.push_back(m_Labels[j]);
		Label2s.push_back(m_Labels[i]);
		Scores.push_back(m_Scores[k]);
		}

	m_SB.m_SBS = SBS_OtherAlgoScore;
	m_SB.SetHits(Label1s, Label2s, Scores);
	m_SB.SetScoreOrder();

	string Msg;
	Ps(Msg, "%s %s %s gap %d/%d\n",
		m_AlignMethod.c_str(),
		m_SubstMxName.c_str(),
		m_ByteSeqMethod.c_str(),
		Paralign::m_Open,
		Paralign::m_Ext);

	m_SB.WriteOutput(Msg);
	}

void ParaSearch::WriteHits(const string &FN) const
	{
	if (FN == "")
		return;

	FILE *f = CreateStdioFile(FN);
	for (uint k = 0; k < m_PairCount; ++k)
		{
		uint i, j;
		triangle_k_to_ij(k, m_SeqCount, i, j);
		fprintf(f, "%.3g", m_Scores[k]);
		fprintf(f, "\t%s", m_Labels[i].c_str());
		fprintf(f, "\t%s", m_Labels[j].c_str());
		fprintf(f, "\n");

		fprintf(f, "%.3g", m_Scores[k]);
		fprintf(f, "\t%s", m_Labels[j].c_str());
		fprintf(f, "\t%s", m_Labels[i].c_str());
		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	}

void ParaSearch::SetDomIdxs()
	{
	asserta(SIZE(m_Labels) == m_SeqCount);
	m_DomIdxs.clear();
	m_DomIdxs.reserve(m_SeqCount);
	for (uint i = 0; i < m_SeqCount; ++i)
		m_DomIdxs.push_back(m_SB.GetDomIdx(m_Labels[i]));
	}

void ParaSearch::ClearHitsAndResults()
	{
	Paralign::ClearStats();
	m_SB.ClearHitsAndResults();
	}

void ParaSearch::SetGapParams(int Open, int Ext)
	{
	Paralign::m_Open = Open;
	Paralign::m_Ext = Ext;
	}

// -seqsmethod		mu | numu (also mux but redundant)
// -alignmethod		para | sw
// -mxname			Mu_S_k_i8 | Mu_scop40_tm0_6_0_8_fa2 | musubstmx
void cmd_para_scop40()
	{
	ParaSearch PS;
	PS.ReadLookup(opt(lookup));
	PS.GetByteSeqs(g_Arg1, opt(seqsmethod));
	PS.SetDomIdxs();
	Paralign::SetSubstMxByName(opt(mxname));
	PS.Search(opt(alignmethod));
	PS.WriteHits(opt(output));
	PS.Bench();
	}