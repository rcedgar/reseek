#include "myutils.h"
#include "scop40bench.h"
#include "seqdb.h"
#include "alpha.h"
#include "cigar.h"
#include "sort.h"
#include "binner.h"

#define ALLOW_MISSING_DOMS	0

static FILE *s_fIn;
static FILE *s_fOut;
static mutex s_InputLock;
static mutex s_OutpuLock;
static uint s_PairCount;
static uint s_NextPairIndex;
static bool s_ShowProgress = true;
static time_t s_timeLastProgress;

void SCOP40Bench::StaticThreadBodyRealign(uint ThreadIndex, SCOP40Bench *ptrSB)
	{
	ptrSB->ThreadBodyRealign(ThreadIndex);
	}

void SCOP40Bench::ThreadBodyRealign(uint ThreadIndex)
	{
	DSSAligner DA;
	DA.SetParams(*m_Params);

	string Line;
	vector<string> Fields;
	for (;;)
		{
		bool Ok = false;
		s_InputLock.lock();
		uint PairIndex = s_NextPairIndex;
		if (s_NextPairIndex < s_PairCount)
			{
			if (s_ShowProgress && PairIndex%1000 == 0)
				{
				time_t now = time(0);
				if (now > s_timeLastProgress)
					{
					Progress("%u / %u\r", PairIndex, s_PairCount);
					s_timeLastProgress = now;
					}
				}
			++s_NextPairIndex;
			}
		s_InputLock.unlock();
		if (PairIndex == s_PairCount)
			return;

		uint ChainIdx1 = m_RealChainIdxs1[PairIndex];
		uint ChainIdx2 = m_RealChainIdxs2[PairIndex];
#if ALLOW_MISSING_DOMS
		if (ChainIdx1 == UINT_MAX || ChainIdx2 == UINT_MAX)
			continue;
#endif

		const PDBChain &Chain1 = *m_DBChains[ChainIdx1];
		const PDBChain &Chain2 = *m_DBChains[ChainIdx2];

		const vector<vector<byte> > &Profile1 = *m_DBProfiles[ChainIdx1];
		const vector<vector<byte> > &Profile2 = *m_DBProfiles[ChainIdx2];

		DA.SetQuery(Chain1, &Profile1, 0, 0, 0);
		DA.SetTarget(Chain2, &Profile2, 0, 0, 0);
		DA.Align_SWOnly();
		m_Scores[PairIndex] = DA.m_AlnFwdScore;
		}
	}

void SCOP40Bench::SetupRealign(const string &TsvFN)
	{
	m_DomToChainIdx.clear();
	const uint ChainCount = GetDBChainCount();
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const string &Label = m_DBChains[ChainIdx]->m_Label;
		string Dom;
		GetDomFromLabel(Label, Dom);
		m_DomToChainIdx[Dom] = ChainIdx;
		}

	FILE *fTsv = OpenStdioFile(TsvFN);
	s_PairCount = 0;
	s_NextPairIndex = 0;

	asserta(m_RealChainIdxs1.empty());
	asserta(m_RealChainIdxs2.empty());
	asserta(m_RealTs.empty());
	m_RealNT = 0;
	m_RealNF = 0;
	string Line;
	vector<string> Fields;
	for (;;)
		{
		bool Ok = ReadLineStdioFile(fTsv, Line);
		if (!Ok)
			break;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 6);
		const string &Dom1 = Fields[0];
		const string &Scopid1 = Fields[1];
		const string &Dom2 = Fields[2];
		const string &Scopid2 = Fields[3];
		float DP = (float) StrToFloat(Fields[4]);
		const string &sTF = Fields[5];

		map<string, uint>::const_iterator iter1 = m_DomToChainIdx.find(Dom1);
		map<string, uint>::const_iterator iter2 = m_DomToChainIdx.find(Dom2);
		asserta(iter1 != m_DomToChainIdx.end());
		asserta(iter2 != m_DomToChainIdx.end());

		uint ChainIdx1 = iter1->second;
		uint ChainIdx2 = iter2->second;

		const string &ChainLabel1 = m_DBChains[ChainIdx1]->m_Label;
		const string &ChainLabel2 = m_DBChains[ChainIdx2]->m_Label;

		string ChainDom1;
		string ChainDom2;
		GetDomFromLabel(ChainLabel1, ChainDom1);
		GetDomFromLabel(ChainLabel2, ChainDom2);
		asserta(ChainDom1 == Dom1);
		asserta(ChainDom2 == Dom2);

		string SF1, SF2;
		GetSFFromScopid(Scopid1, SF1);
		GetSFFromScopid(Scopid2, SF2);

		bool T = (SF1 == SF2);
		m_RealTs.push_back(T);
		if (T == 1)
			++m_RealNT;
		else
			++m_RealNF;
		m_RealChainIdxs1.push_back(ChainIdx1);
		m_RealChainIdxs2.push_back(ChainIdx2);
		m_RealTs.push_back(T);
		}
	s_PairCount = SIZE(m_RealChainIdxs1);
	m_Scores.resize(s_PairCount, FLT_MAX);
	}

void SCOP40Bench::RunRealign()
	{
	const uint ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	time_t t0 = time(0);
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadBodyRealign, ThreadIndex, this);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	asserta(SIZE(m_Scores) == s_PairCount);
	time_t t1 = time(0);
	AnalyzeRealign();
	}

void SCOP40Bench::AnalyzeRealign()
	{
	asserta(SIZE(m_Scores) == s_PairCount);
	vector<uint> Order(UINT_MAX);
	QuickSortOrder(m_Scores.data(), s_PairCount, Order.data());
	uint nt = 0;
	uint nf = 0;
	uint best_ne = UINT_MAX;
	uint best_tp = UINT_MAX;
	uint best_fp = UINT_MAX;
	uint best_tn = UINT_MAX;
	uint best_fn = UINT_MAX;
	float BestScore = FLT_MAX;
	float MaxFPScore = -FLT_MAX;
	float MaxScore = -FLT_MAX;
	vector<float> TPScores;
	vector<float> FPScores;
	float LastScore = -FLT_MAX;
	for (uint k = 0; k < s_PairCount; ++k)
		{
		uint PairIndex = Order[k];

		float Score = m_Scores[PairIndex];
		uint T = m_RealTs[PairIndex];

		asserta(Score >= LastScore);
		LastScore = Score;
		MaxScore = max(Score, MaxScore);
		if (T == 0)
			{
			TPScores.push_back(Score);
			++nf;
			}
		else if (T == 1)
			{
			FPScores.push_back(Score);
			MaxFPScore = max(Score, MaxFPScore);
			++nt;
			}
		else
			asserta(false);
		asserta(nt <= m_RealNT);
		asserta(nf <= m_RealNF);
		uint nfn = m_RealNT - nt;
		uint ne = nfn + nf;
		if (ne < best_ne)
			{
			best_ne = ne;
			best_tp = nt;
			best_fp = nf;
			best_fn = nfn;
			best_tn = m_RealNF - nf;
			BestScore = Score;
			}
		}
	ProgressLog("ne %u, tp %u, fp %u, fn %u, tn %u, score %.3g\n",
				best_ne,
				best_tp,
				best_fp,
				best_fn,
				best_tn,
				BestScore);

	float MaxBinScore = min(MaxScore, 2*MaxFPScore);
	Binner<float> TPB(TPScores, 32, 0, MaxBinScore);
	Binner<float> FPB(FPScores, 32, 0, MaxBinScore);
	for (uint Bin = 0; Bin < 32; ++Bin)
		{
		uint TPCount = TPB.GetCount(Bin);
		uint FPCount = FPB.GetCount(Bin);
		float Mid = TPB.GetBinMid(Bin);
		ProgressLog("%10u  %10u  %.3g\n", TPCount, FPCount, Mid);
		}
	}

void cmd_scop40bench_real()
	{
	optset_fast = true;
	opt_fast = true;

	const string CalFN = g_Arg1;
	s_fOut = CreateStdioFile(opt(output));

	DSSParams Params;
	Params.SetDSSParams(DM_UseCommandLineOption, SCOP40_DBSIZE);

	SCOP40Bench SB;
	SB.m_Params = &Params;
	SB.LoadDB(CalFN);
	SB.SetupRealign(opt(input2));
	
	float MaxFPR = 0.005f;
	if (optset_maxfpr)
		MaxFPR = (float) opt(maxfpr);

	SB.m_QuerySelf = true;
	SB.m_ScoresAreEvalues = true;
	if (opt(scores_are_not_evalues))
		SB.m_ScoresAreEvalues = false;
	SB.RunRealign();
	CloseStdioFile(s_fOut);
	s_fOut = 0;
	}
