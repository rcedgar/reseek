#include "myutils.h"
#include "parasearch.h"
#include "triangle.h"
#include "nu.h"

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

void ParaSearch::Search(const string &AlignMethod, string SubstMxName)
	{
	m_AlignMethod = AlignMethod;
	m_SubstMxName = SubstMxName;

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
	Paralign::SetSubstMx(SubstMxName);
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
	// Make sure ThreadIdx is defined inside the parallel region if used locally
	uint ThreadIdx = GetThreadIndex(); // Or your GetThreadIndex()
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
	// Threads implicitly wait at the end of the omp for construct
	}

	ProgressStep(m_PairCount-1, m_PairCount, "Aligning");
	ProgressLog("%u long, %u saturated, %u 8-bit, %u 16-bit, %u SW\n",
		Paralign::m_TooLongCount.load(),
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

void ParaSearch::GetByteSeqs_numu(const string &FN)
	{
	Nu A;
	A.SetMu();

	ReadChains(FN, m_Chains);
	const uint ChainCount = SIZE(m_Chains);
	vector<vector<byte> > ByteSeqs(ChainCount);
	vector<string> Labels;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Nu::GetLetters(Mu)");
		const PDBChain &Chain = *m_Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		vector<byte> &ByteSeq = ByteSeqs[ChainIdx];
		A.GetLetters(Chain, ByteSeq);
		if (opt(fixmubyteseq))
			FixMuByteSeq(ByteSeq);
		Labels.push_back(Chain.m_Label);
		}
	}

void ParaSearch::Bench(const string &LookupFN)
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
		}
	m_SB.SetHits(Label1s, Label2s, Scores);
	m_SB.m_SBS = SBS_OtherAlgoScore;
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

	Die("TODO");
	//const uint HitCount = SIZE(m_Label1s);
	//FILE *f = CreateStdioFile(FN);
	//for (uint i = 0; i < HitCount; ++i)
	//	{
	//	ProgressStep(i, HitCount, "Writing %s", FN.c_str());

	//	fprintf(f, "%s\t%s\t%.3g\n",
	//		m_Label1s[i].c_str(),
	//		m_Label2s[i].c_str(),
	//		m_Scores[i]);
	//	}
	//CloseStdioFile(f);
	}

void ParaSearch::SetDomIdxs()
	{
	asserta(SIZE(m_Labels) == m_SeqCount);
	m_DomIdxs.clear();
	m_DomIdxs.reserve(m_SeqCount);
	for (uint i = 0; i < m_SeqCount; ++i)
		m_DomIdxs.push_back(m_SB.GetDomIdx(m_Labels[i]));
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
	PS.Search(opt(alignmethod), opt(mxname));
	PS.WriteHits(opt(output));
	PS.Bench(opt(lookup));
	}