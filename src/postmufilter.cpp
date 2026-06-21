#if 0
/***
Query bags are pre-computed.
ThreadBody()
	1. Load next target in DB order.
	2. Create target bag.
	3. For each query with prefilter hit:
		3b. Align to query
***/

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

static atomic<uint> s_TargetCounter;
static uint s_TargetCount;
static const vector<ChainBag *> *s_ptrCBQs;
static const vector<uint> *s_ptrTargetIdxs;
static const unordered_map<uint, vector<uint> > *s_ptrTargetIdxToQueryIdxs;
static BCAData *s_ptrDB;
static uint s_ScannedCount;
static double s_MaxEvalue = 10;
static double s_MaxPvalue = -1;
static double s_MinTS = 9e9;
static FILE *s_fTsv;
static FILE *s_fAln;

static bool Accept(const DSSAligner &DA)//TODO compare Search
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
		ProgressStep(0, s_TargetCount, "Scanning");

	for (;;)
		{
		uint TargetCounter = s_TargetCounter.fetch_add(1, std::memory_order_relaxed);
		if (TargetCounter >= s_TargetCount)
			{
			if (ThreadIndex == 0)
				ProgressStep(s_TargetCount-1, s_TargetCount, "Scanning");
			return;
			}
		if (ThreadIndex == 0 && TargetCounter + 1 < s_TargetCount)
			ProgressStep(TargetCounter, s_TargetCount, "Scanning");

		const uint TargetIdx = (*s_ptrTargetIdxs)[TargetCounter];

		PDBChain ChainT;
		DB.ReadChain(TargetIdx, ChainT);
		MakeBag(ChainT, CBT, D, DASelfRevT, MKF);

		unordered_map<uint, vector<uint> >::const_iterator iter =
			s_ptrTargetIdxToQueryIdxs->find(TargetIdx);
		asserta(iter != s_ptrTargetIdxToQueryIdxs->end());
		const vector<uint> &QueryIdxs = iter->second;

		const uint QN = SIZE(QueryIdxs);
		for (uint k = 0; k < QN; ++k)
			{
			if (ThreadIndex == 0 && TargetCounter + 1 < s_TargetCount)
				ProgressStep(TargetCounter, s_TargetCount, "Scanning");

			const uint QueryIdx = QueryIdxs[k];
			const ChainBag &CBQ = *CBQs[QueryIdx];
			const PDBChain &ChainQ = *CBQ.m_ptrChain;
			TheDA.m_MKF.SetQ(ChainQ.m_Label, CBQ.m_ptrMuLetters, CBQ.m_ptrMuKmers);
			TheDA.SetBagA(CBQ);
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
void PostMuFilter(
	const vector<ChainBag *> &CBQs,
	const string &DBBCAFN,
	const vector<uint> &TargetIdxs,
	const unordered_map<uint, vector<uint> > &TargetIdxToQueryIdxs,
	const string &HitsFN)
	{
	time_t t0 = time(0);

	s_ptrCBQs = &CBQs;
	s_ptrTargetIdxs = &TargetIdxs;
	s_ptrTargetIdxToQueryIdxs = &TargetIdxToQueryIdxs;

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
	
	s_TargetCount = SIZE(TargetIdxs);
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
	ProgressLog("%10u  m_PostMuFilterMKFCount\n", DSSAligner::m_PostMuFilterMKFCount.load());
	ProgressLog("%10u  m_PostMuFilterOmegaDiscardCount\n", DSSAligner::m_PostMuFilterOmegaDiscardCount.load());
	ProgressLog("%10u  m_PostMuFilterSWCount\n", DSSAligner::m_PostMuFilterSWCount.load());
	ProgressLog("%10u  m_XDropDiscardCount1\n", DSSAligner::m_XDropDiscardCount1.load());
	ProgressLog("%10u  m_XDropDiscardCount2\n", DSSAligner::m_XDropDiscardCount2.load());
	ProgressLog("%10u  m_XDropAlnCount\n", DSSAligner::m_XDropAlnCount.load());
	}
#endif // 0