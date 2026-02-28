/***
Query bags are pre-computed.
ThreadBody()
	1. Load next query
	2. For each target with prefilter hit:
		2a. Create target bag
		2b. Align to query
Much slower for all-vs-all e.g. SCOP40
***/
#if 0
#include "myutils.h"
#include "dssaligner.h"
#include "chainreader2.h"
#include "alncounts.h"
#include "prefiltermu.h"

float GetSelfRevScore(DSSAligner &DA, DSS &D, const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers);
void MakeBags(const vector<PDBChain *> Chains, vector<ChainBag *> &Bags);
void MakeBag(const PDBChain &Chain, ChainBag &CB,
	DSS &D, DSSAligner &DASelfRev, MuKmerFilter &MKF);

static atomic<uint> s_QueryIdx;
static uint s_QueryCount;
static const vector<ChainBag *> *s_ptrCBQs;
static BCAData *s_ptrDB;
static uint s_ScannedCount;
static double s_MaxEvalue = 10;
static double s_MaxPvalue = -1;
static double s_MinTS = 9e9;
static FILE *s_fTsv;
static FILE *s_fAln;

static bool Accept(const DSSAligner &DA)//@@TODO compare Search
	{
	if (DA.m_EvalueA <= s_MaxEvalue)
		return true;
	if (DA.m_PvalueA <= s_MaxPvalue)
		return true;
	if (DA.m_NewTestStatisticA >= s_MinTS)
		return true;
	return false;
	}

static void ThreadBody_Scan(uint ThreadIndex)
	{
	DSS D;
	DSSAligner DASelfRevQ;
	DSSAligner DASelfRevT;
	MuKmerFilter MKF;
	const vector<ChainBag *> &CBQs = *s_ptrCBQs;
	const BCAData &DB = *s_ptrDB;
	vector<vector<byte> > DBProfile;
	vector<byte> DBMuLetters;
	vector<uint> DBMuKmers;
	float SelfRevScore = 0;
	ChainBag CBT;
	DSSAligner TheDA;

	string Line;
	vector<string> Fields;

	if (ThreadIndex == 0)
		ProgressStep(0, s_QueryCount, "Scanning");

	for (;;)
		{
		uint QueryIdx = s_QueryIdx.fetch_add(1, std::memory_order_relaxed);
		if (QueryIdx >= s_QueryCount)
			{
			if (ThreadIndex == 0)
				ProgressStep(s_QueryCount-1, s_QueryCount, "Scanning");
			return;
			}
		if (ThreadIndex == 0 && QueryIdx + 1 < s_QueryCount)
			ProgressStep(QueryIdx, s_QueryCount, "Scanning");

		const ChainBag &CBQ = *CBQs[QueryIdx];
		const PDBChain &ChainQ = *CBQ.m_ptrChain;
		TheDA.m_MKF.SetQ(ChainQ.m_Label, CBQ.m_ptrMuLetters, CBQ.m_ptrMuKmers);
		TheDA.SetBagA(CBQ);

		const vector<uint> &TargetIdxs = PrefilterMu::m_RSB.GetTargetIdxs(QueryIdx);
		const uint TIN = SIZE(TargetIdxs);
		for (uint k = 0; k < TIN; ++k)
			{
			if (ThreadIndex == 0 && QueryIdx + 1 < s_QueryCount)
				ProgressStep(QueryIdx, s_QueryCount, "Scanning");

			uint TargetIdx = TargetIdxs[k];

			PDBChain DBChain;
			DB.ReadChain(TargetIdx, DBChain);

			D.Init(DBChain);
			D.GetProfile(DBProfile);
			D.GetMuLetters(DBMuLetters);
			D.GetMuKmers(DBMuLetters, DBMuKmers, DSSParams::m_MKFPatternStr);

			//@@TODO usually don't need full self rev score
			float DBSelfRevScore = GetSelfRevScore(DASelfRevT, D, DBChain, DBProfile,
											   &DBMuLetters, &DBMuKmers);

			CBT.m_ptrChain = &DBChain;
			CBT.m_ptrProfile = &DBProfile;
			CBT.m_ptrMuLetters = &DBMuLetters;
			CBT.m_ptrMuKmers = &DBMuKmers;
			CBT.m_SelfRevScore = DBSelfRevScore;
			CBT.m_ptrProfPara8 = DASelfRevT.m_ProfPara8;
			CBT.m_ptrProfPara16 = DASelfRevT.m_ProfPara16;
			CBT.m_ptrProfParaRev8 = DASelfRevT.m_ProfParaRev8;
			CBT.m_ptrProfParaRev16 = DASelfRevT.m_ProfParaRev16;
			TheDA.AlignBagB(CBT);
			if (Accept(TheDA))
				{
				incac(scanhits);
				TheDA.ToTsv(s_fTsv, true);
				TheDA.ToAln(s_fAln, true);
				}
			else
				incac(scanrejects);
			DASelfRevT.UnsetQuery();
			}
		}
	}

// Query & DB need C-alpha
void PostMuFilter(const vector<ChainBag *> &CBQs,
				  const string &DBBCAFN,
				  const string &HitsFN)
	{
	time_t t0 = time(0);
	s_ptrCBQs = &CBQs;
	if (optset_evalue)
		s_MaxEvalue = opt(evalue);
	else if (optset_verysensitive)
		s_MaxEvalue = 9e9;
	if (optset_pvalue)
		s_MaxPvalue = opt(pvalue);
	if (optset_mints)
		s_MinTS = opt(mints);

	s_fAln = CreateStdioFile(opt(aln));
	s_fTsv = CreateStdioFile(HitsFN);
	
	//vector<PDBChain *> QChains;
	//ReadChains(QueryCAFN, QChains);
	s_QueryCount = SIZE(CBQs);
	setac(queries, s_QueryCount);

	BCAData DB;
	DB.Open(DBBCAFN);
	s_ptrDB = &DB;
	setac(targets, DB.GetChainCount());

	uint ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;

	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody_Scan, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	CloseStdioFile(s_fAln);
	CloseStdioFile(s_fTsv);
	time_t t1 = time(0);
	ProgressLog("Post-mu %u secs\n", uint(t1 - t0));
	ProgressLog("%10u  m_AlignBagB_MKFCount\n", DSSAligner::m_AlignBagB_MKFCount.load());
	ProgressLog("%10u  m_PostMuFilterOmegaDiscardCount\n", DSSAligner::m_PostMuFilterOmegaDiscardCount.load());
	ProgressLog("%10u  m_PostMuFilterSWCount\n", DSSAligner::m_PostMuFilterSWCount.load());
	ProgressLog("%10u  m_XDropDiscardCount1\n", DSSAligner::m_XDropDiscardCount1.load());
	ProgressLog("%10u  m_XDropDiscardCount2\n", DSSAligner::m_XDropDiscardCount2.load());
	ProgressLog("%10u  m_XDropAlnCount\n", DSSAligner::m_XDropAlnCount.load());
	}
#endif