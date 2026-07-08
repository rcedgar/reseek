#if 0
#include "myutils.h"
#include "dss.h"
#include "seqdb.h"
#include "museqsource.h"
#include "dbsearcher.h"
#include "output.h"
#include "statsig.h"
#include "alpha.h"

void MakeBags(const vector<PDBChain *> Chains, vector<ChainBag *> &Bags);
void MuPreFilter(SeqDB &QDB, MuSeqSource &FSS, vector<uint> &TargetIdxs,
	unordered_map<uint, vector<uint> > &TargetIdxToQueryIdxs);

void PostMuFilter(
	const vector<ChainBag *> &CBQs,
	const string &DBBCAFN,
	const vector<uint> &TargetIdxs,
	const unordered_map<uint, vector<uint> > &TargetIdxToQueryIdxs,
	const string &HitsFN);

void SelfSearch()
	{
	const string &QFN = g_Arg1;
	if (optset_db)
		Die("-db not used for -selfsearch");

	DBSearcher DBS;

	DBS.LoadDB(QFN);
	DBS.Setup();

	OpenOutputFiles();
	DBS.RunSelf();
	CloseOutputFiles();
	}

static void Search_NoMuFilter()
	{
	if (!optset_db)
		Die("-db required");

	const string &QFN = g_Arg1;
	const string &DBFN = opt(db);

	DBSearcher DBS;

	DBS.LoadDB(QFN);
	DBS.Setup();

	OpenOutputFiles();
	ChainReader2 CR;
	CR.Open(DBFN);
	DBS.RunQuery(CR);
	CloseOutputFiles();
	}

void MakeMuSeqDB(const vector<ChainBag *> &CBs, SeqDB &DB)
	{
	DB.Clear();
	const uint ChainCount = SIZE(CBs);
	DB.m_Labels.reserve(ChainCount);
	DB.m_Seqs.resize(ChainCount);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const ChainBag &CB = *CBs[ChainIdx];
		const vector<byte> &MuLetters = *CB.m_ptrMuLetters;
		const uint L = CB.GetSeqLength();
		DB.m_Labels.push_back(CB.GetLabel());
		string &Seq = DB.m_Seqs[ChainIdx];
		Seq.reserve(L);
		for (uint i = 0; i < L; ++i)
			{
			assert(MuLetters[i] < 36);
			Seq += g_LetterToCharMu[MuLetters[i]];
			}
		}
	}

void cmd_mufilter()
	{
	const string &QFN = g_Arg1;

	DSSParams::Init(DM_UseCommandLineOption);
	DBSearcher DBS;

	DBS.LoadDB(QFN);
	DBS.Setup();

	OpenOutputFiles();
	DBS.RunSelf();
	CloseOutputFiles();
	}

void cmd_search()
	{
	DSSParams::Init(DM_UseCommandLineOption);
	if (!optset_db)
		{
		SelfSearch();
		return;
		}

	if (!optset_fast)
		{
		Search_NoMuFilter();
		return;
		}

	const string &QueryFN = g_Arg1;
	const string &DBFN = string(opt(db));

	if (!EndsWith(DBFN, ".bca"))
		Die(".bca format required for -db");

	vector<PDBChain *> ChainsQ;
	ReadChains(QueryFN, ChainsQ);

	vector<ChainBag *> CBQs;
	MakeBags(ChainsQ, CBQs);

	SeqDB MuQueryDB;
	//MuQueryDB.FromSS(QSS);
	MakeMuSeqDB(CBQs, MuQueryDB);

	MuSeqSource DBSS;
	DBSS.OpenChains(DBFN);

	vector<uint> TargetIdxs;
	unordered_map<uint, vector<uint> > TargetIdxToQueryIdxs;
	MuPreFilter(MuQueryDB, DBSS, TargetIdxs, TargetIdxToQueryIdxs);

	DSSParams::SetAlgoMode(DM_AlwaysFast);
	PostMuFilter(CBQs, DBFN, TargetIdxs, TargetIdxToQueryIdxs, opt(output));
	}
#endif
