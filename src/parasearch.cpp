#include "myutils.h"
#include "parasearch.h"
#include "triangle.h"
#include "nu.h"
#include "sort.h"
#include <numeric>

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

void ParaSearch::AppendHit_rev(uint i, uint j, float Score)
	{
	uint k = triangle_ij_to_k(i, j, m_SeqCount);
	m_Scores_rev[k] = Score;
	}

void ParaSearch::AppendHit(uint i, uint j, float Score)
	{
	uint k = triangle_ij_to_k(i, j, m_SeqCount);
	if (m_DoReverse)
		m_Scores_fwd[k] = Score;
	m_Scores[k] = Score;
	}

float ParaSearch::GetSelfScore_rev(Paralign &PA, uint ChainIdx)
	{
	asserta(m_AlignMethod == "para");
	const string &Label = m_Labels[ChainIdx];
	const vector<byte> &ByteSeq = m_ByteSeqs[ChainIdx];
	uint L = SIZE(ByteSeq);
	PA.m_Q = ByteSeq.data();
	PA.m_LQ = L;

	PA.SetQueryProfile(Label, ByteSeq.data(), L);
	PA.Align_ScoreOnly_rev(Label, ByteSeq.data(), L);
	return (float) PA.m_Score_rev;
	}

void ParaSearch::SetSelfScores_rev(const string &AlignMethod)
	{
	InitThreads(AlignMethod, true);
	if (m_SelfScores_rev != 0)
		myfree(m_SelfScores_rev);
	m_SelfScores_rev = myalloc(float, m_SeqCount);
	Paralign PA;
	PA.m_DoReverse = true;
	Progress("Self scores... ");
	for (uint i = 0; i < m_SeqCount; ++i)
		m_SelfScores_rev[i] = GetSelfScore_rev(PA, i);
	Progress("done\n");
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
		asserta(!m_DoReverse);
		PA.Align_SWFast(Label_j, ByteSeq_j.data(), L_j);
		AppendHit(i, j, PA.m_SWFastScore);
		}
	else if (m_AlignMethod == "para")
		{
		PA.Align_ScoreOnly(Label_j, ByteSeq_j.data(), L_j);
		AppendHit(i, j, (float) PA.m_Score);
		if (m_DoReverse)
			{
			PA.Align_ScoreOnly_rev(Label_j, ByteSeq_j.data(), L_j);
			AppendHit_rev(i, j, (float) PA.m_Score_rev);
			}
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
		{
		PA.SetQueryProfile(Label_i, ByteSeq_i.data(), L_i);
		if (m_DoReverse)
			PA.SetQueryProfile_rev(ByteSeq_i.data(), L_i);
		}
	else if (m_AlignMethod == "sw")
		;
	else
		Die("m_AlignMethod=%s", m_AlignMethod.c_str());
	}

void ParaSearch::InitThreads(const string &AlignMethod, bool DoReverse)
	{
	m_AlignMethod = AlignMethod;

	ProgressLog("Search %s %s %s\n",
		m_AlignMethod.c_str(),
		m_SubstMxName.c_str(),
		m_ByteSeqMethod.c_str());

	asserta(m_SeqCount == SIZE(m_Chains));
	asserta(m_SeqCount == SIZE(m_ByteSeqs));
	uint PairCount2 = triangle_get_K(m_SeqCount) + 1;
	asserta(m_PairCount == PairCount2);
	if (m_Scores != 0)
		myfree(m_Scores);
	if (m_Scores_fwd != 0)
		myfree(m_Scores_fwd);
	if (m_Scores_rev != 0)
		myfree(m_Scores_rev);
	m_Scores = myalloc(float, m_PairCount);
	if (m_DoReverse)
		{
		m_Scores_fwd = myalloc(float, m_PairCount);
		m_Scores_rev = myalloc(float, m_PairCount);
		}
	const uint ThreadCount = GetRequestedThreadCount();
	m_PAs.clear();
	m_QueryIdxs.clear();
	for (uint i = 0; i < ThreadCount; ++i)
		{
		Paralign *PA = new Paralign;
		PA->m_DoReverse = DoReverse;
		m_PAs.push_back(PA);
		m_QueryIdxs.push_back(UINT_MAX);
		}
	}

void ParaSearch::Search(const string &AlignMethod, bool DoReverse)
	{
	InitThreads(AlignMethod, DoReverse);
	atomic<uint> Counter = 0;
	ProgressStep(0, m_PairCount, "Aligning");

	const uint ThreadCount = GetRequestedThreadCount();
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

void ParaSearch::SetLookupFromLabels()
	{
	m_SFIdxToSize.clear();
	m_SFIdxToSize.resize(2000, 0);

	const uint N = SIZE(m_Labels);
	for (uint LabelIdx = 0; LabelIdx < N; ++LabelIdx)
		{
		const string &Label = m_Labels[LabelIdx];
		vector<string> Fields, Fields2;
		Split(Label, Fields, '/');
		asserta(SIZE(Fields) == 2);
		const string &Dom = Fields[0];
		const string &ScopId = Fields[1];
		Split(ScopId, Fields2, '.');
		asserta(SIZE(Fields2) == 4);
		const string &SF = Fields2[0] + "." + Fields2[1] + "." + Fields2[2];
		AddDom(Dom, SF, LabelIdx);
		}

	m_NT = 0;
	m_NF = 0;
	const uint SFCount = SIZE(m_SFs);
	for (uint SFIdx = 0; SFIdx < SFCount; ++SFIdx)
		{
		uint Size = m_SFIdxToSize[SFIdx];
		asserta(Size > 0);
		m_NT += (Size*(Size - 1))/2;
		}
	m_NT *= 2;
	uint DomCount = SIZE(m_Labels);
	m_NF = DomCount*(DomCount-1) - m_NT;
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

void ParaSearch::SetScoreOrder()
	{
	uint K = triangle_get_K(m_SeqCount);
	if (m_ScoreOrder == 0)
		myfree(m_ScoreOrder);
	m_ScoreOrder = myalloc(uint, K);
	QuickSortOrderDesc(m_Scores, K, m_ScoreOrder);
	}

void ParaSearch::BenchRev(const string &Msg, 
	float SelfWeight, float RevWeight)
	{
	asserta(m_DoReverse);
	asserta(m_SelfScores_rev != 0);
	asserta(m_Scores_fwd != 0);
	asserta(m_Scores_rev != 0);
	asserta(m_Scores != 0);

	uint K = triangle_get_K(m_SeqCount);
	for (uint HitIdx = 0; HitIdx < K; ++HitIdx)
		{
		uint LabelIdx_i, LabelIdx_j;
		triangle_k_to_ij(HitIdx, m_SeqCount, LabelIdx_i, LabelIdx_j);
		if (LabelIdx_i == LabelIdx_j)
			continue;
		float SelfScore_rev_i = m_SelfScores_rev[LabelIdx_i];
		float SelfScore_rev_j = m_SelfScores_rev[LabelIdx_j];
		float Score_fwd = m_Scores_fwd[HitIdx];
		float Score_rev = m_Scores_rev[HitIdx];
		
		m_Scores[HitIdx] = Score_fwd - RevWeight*Score_rev -
			SelfWeight*(SelfScore_rev_i + SelfScore_rev_j);
		}
	SetScoreOrder();
	Bench(Msg);
	}

void ParaSearch::Bench(const string &Msg)
	{
	asserta(m_ScoreOrder != 0);
	uint K = triangle_get_K(m_SeqCount);
	uint nt = 0;
	uint nf = 0;
	float LastScore = FLT_MAX;
	float SEPQ0_1 = FLT_MAX;
	float SEPQ1 = FLT_MAX;
	float SEPQ10 = FLT_MAX;
	for (uint k = 0; k < K; ++k)
		{
		uint HitIdx = m_ScoreOrder[k];
		uint LabelIdx_i, LabelIdx_j;
		triangle_k_to_ij(HitIdx, m_SeqCount, LabelIdx_i, LabelIdx_j);
		if (LabelIdx_i == LabelIdx_j)
			continue;
		float Score = m_Scores[HitIdx];
		if (Score != LastScore)
			{
			asserta(Score < LastScore);
			float EPQ = 2*float(nf)/m_SeqCount;
			float Sens = 2*float(nt)/m_NT;
			if (SEPQ0_1 == FLT_MAX && EPQ >= 0.1) SEPQ0_1 = Sens;
			if (SEPQ1 == FLT_MAX   && EPQ >= 1)   SEPQ1   = Sens;
			if (SEPQ10 == FLT_MAX  && EPQ >= 10)  SEPQ10  = Sens;
			LastScore = Score;
			}
		uint SFIdx_i = m_LabelIdxToSFIdx[LabelIdx_i];
		uint SFIdx_j = m_LabelIdxToSFIdx[LabelIdx_j];
		if (SFIdx_i == SFIdx_j)
			++nt;
		else
			++nf;
		}
	float EPQ = float(nf)/m_SeqCount;
	float Sens = float(nt)/m_NT;
	if (SEPQ0_1 == FLT_MAX && EPQ >= 0.1) SEPQ0_1 = Sens;
	if (SEPQ1 == FLT_MAX   && EPQ >= 1)   SEPQ1   = Sens;
	if (SEPQ10 == FLT_MAX  && EPQ >= 10)  SEPQ10  = Sens;
	m_Sum3 = SEPQ0_1*2 + SEPQ1*3/2 + SEPQ10;

	if (Msg != "")
		ProgressLog("%s ", Msg.c_str());
	ProgressLog("%s %s %s gap %d/%d N=%u NT=%u\n",
		m_AlignMethod.c_str(),
		m_SubstMxName.c_str(),
		m_ByteSeqMethod.c_str(),
		Paralign::m_Open,
		Paralign::m_Ext,
		m_SeqCount,
		m_NT);
	ProgressLog("SEPQ0.1=%.3f", SEPQ0_1);
	ProgressLog(" SEPQ1=%.3f", SEPQ1);
	ProgressLog(" SEPQ10=%.3f", SEPQ10);
	ProgressLog(" Sum3=%.3f", m_Sum3);
	ProgressLog("\n");
	}

void ParaSearch::AddDom(
	const string &Dom, const string &SF, uint LabelIdx)
	{
	uint SFIdx = UINT_MAX;
	if (m_SFToIdx.find(SF) == m_SFToIdx.end())
		{
		SFIdx = SIZE(m_SFs);
		m_SFs.push_back(SF);
		m_SFToIdx[SF] = SFIdx;
		}
	else
		SFIdx = m_SFToIdx[SF];
 
	m_LabelIdxToSFIdx.push_back(SFIdx);

	asserta(SFIdx < SIZE(m_SFIdxToSize));
	m_SFIdxToSize[SFIdx] += 1;
	}

void ParaSearch::WriteRevTsv(const string &FN) const
	{
	asserta(m_DoReverse);
	if (FN == "")
		return;
	asserta(m_ScoreOrder != 0);
	asserta(m_SelfScores_rev != 0);
	asserta(m_Scores_rev != 0);
	asserta(m_Scores_fwd != 0);

	FILE *f = CreateStdioFile(FN);

	fprintf(f, "%u\n", m_SeqCount);
	for (uint i = 0; i < m_SeqCount; ++i)
		fprintf(f, "%u\t%s\t%.3g\n",
			i, m_Labels[i].c_str(), m_SelfScores_rev[i]);

	uint K = triangle_get_K(m_SeqCount);
	for (uint k = 0; k < K; ++k)
		{
		ProgressStep(k, K, "Writing %s", FN.c_str());
		uint HitIdx = m_ScoreOrder[k];
		uint i, j;
		triangle_k_to_ij(HitIdx, m_SeqCount, i, j);
		if (i == j)
			continue;

		fprintf(f, "%.3g", m_Scores_fwd[HitIdx]);
		fprintf(f, "\t%.3g", m_Scores_rev[HitIdx]);
		fprintf(f, "\t%s", m_Labels[i].c_str());
		fprintf(f, "\t%s", m_Labels[j].c_str());
		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	}

void ParaSearch::WriteHits(const string &FN) const
	{
	if (FN == "")
		return;
	asserta(m_ScoreOrder != 0);

	FILE *f = CreateStdioFile(FN);
	uint K = triangle_get_K(m_SeqCount);
	for (uint k = 0; k < K; ++k)
		{
		ProgressStep(k, K, "Writing %s", FN.c_str());
		uint HitIdx = m_ScoreOrder[k];
		uint i, j;
		triangle_k_to_ij(HitIdx, m_SeqCount, i, j);
		if (i == j)
			continue;

		fprintf(f, "%.3g", m_Scores[HitIdx]);
		fprintf(f, "\t%s", m_Labels[i].c_str());
		fprintf(f, "\t%s", m_Labels[j].c_str());
		fprintf(f, "\n");

		fprintf(f, "%.3g", m_Scores[HitIdx]);
		fprintf(f, "\t%s", m_Labels[j].c_str());
		fprintf(f, "\t%s", m_Labels[i].c_str());
		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	}

void ParaSearch::ClearHitsAndResults()
	{
	Paralign::ClearStats();
	myfree(m_Scores);
	myfree(m_Scores_fwd);
	myfree(m_Scores_rev);
	myfree(m_SelfScores_rev);
	m_Scores = 0;
	m_Scores_fwd = 0;
	m_Scores_rev = 0;
	m_SelfScores_rev = 0;
	m_Sum3 = FLT_MAX;
	for (uint i = 0; i < SIZE(m_PAs); ++i)
		delete m_PAs[i];
	m_PAs.clear();
	}

void ParaSearch::SetGapParams(int Open, int Ext)
	{
	Paralign::m_Open = Open;
	Paralign::m_Ext = Ext;
	}

void ParaSearch::MakeSubset(ParaSearch &Subset, uint SubsetPct)
	{
	vector<uint> ChainIdxs;
	const uint ChainCount = SIZE(m_Labels);
	const uint SubsetChainCount = (ChainCount*SubsetPct)/100;
	asserta(SubsetChainCount > 0 && SubsetChainCount <= ChainCount);
	for (uint i = 0; i < ChainCount; ++i)
		ChainIdxs.push_back(i);
	Shuffle(ChainIdxs);

	Subset.m_SeqCount = SubsetChainCount;
	Subset.m_PairCount = SubsetChainCount*(SubsetChainCount-1)/2 + SubsetChainCount;
	Subset.m_AlignMethod = m_AlignMethod;
	Subset.m_SubstMxName = m_SubstMxName;
	Subset.m_ByteSeqMethod = m_ByteSeqMethod;
	Subset.m_LabelIdxToSFIdx.clear();
	Subset.m_DomIdxs.clear();
	Subset.m_SFs.clear();
	Subset.m_SFIdxToSize.clear();
	Subset.m_NT = UINT_MAX;
	Subset.m_NF = UINT_MAX;
	Subset.m_ScoreOrder = 0;

	Subset.m_Chains.clear();
	Subset.m_Labels.clear();
	Subset.m_ByteSeqs.clear();
	Subset.m_LabelIdxToSFIdx.clear();

	Subset.m_Chains.reserve(SubsetChainCount);
	Subset.m_Labels.reserve(SubsetChainCount);
	Subset.m_ByteSeqs.reserve(SubsetChainCount);
	Subset.m_LabelIdxToSFIdx.reserve(SubsetChainCount);

	for (uint i = 0; i < SubsetChainCount; ++i)
		{
		uint Idx = ChainIdxs[i];
		PDBChain *Chain = m_Chains[Idx];
		Subset.m_Chains.push_back(Chain);
		Subset.m_Labels.push_back(Chain->m_Label);
		Subset.m_ByteSeqs.push_back(m_ByteSeqs[Idx]);
		Subset.m_LabelIdxToSFIdx.push_back(m_LabelIdxToSFIdx[Idx]);
		}
	Subset.SetLookupFromLabels();
	}

// -seqsmethod		mu | numu (also mux but redundant)
// -alignmethod		para | sw
// -mxname			Mu_S_k_i8 | Mu_scop40_tm0_6_0_8_fa2 | musubstmx
void cmd_para_scop40()
	{
	ParaSearch PS;
	PS.GetByteSeqs(g_Arg1, opt(seqsmethod));
	PS.SetLookupFromLabels();
	Paralign::SetSubstMxByName(opt(mxname));
	PS.Search(opt(alignmethod), false);
	PS.WriteHits(opt(output));
	PS.Bench();
	}