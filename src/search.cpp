#include "myutils.h"
#include "dss.h"
#include "seqdb.h"
#include "museqsource.h"
#include "dbsearcher.h"
#include "output.h"
#include "statsig.h"
#include "prefiltermuparams.h"

uint MuPreFilter(const DSSParams &Params,
			  SeqDB &QueryDB,
			  MuSeqSource &FSS,
			  const string &OutputFN);

void PostMuFilter(const DSSParams &Params,
				  const string &MuFilterTsvFN,
				  const string &QueryCAFN,
				  const string &DBBCAFN,
				  const string &HitsFN);

void SelfSearch()
	{
	const string &QFN = g_Arg1;
	if (optset_db)
		Die("-db not used for -selfsearch");

	DBSearcher DBS;
	DSSParams Params;
	Params.SetDSSParams(DM_UseCommandLineOption);
	DBS.m_Params = &Params;

	DBS.LoadDB(QFN);
	DBS.Setup();
	StatSig::Init(DBS.GetDBChainCount());

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
	DSSParams Params;
	Params.SetDSSParams(DM_UseCommandLineOption);
	DBS.m_Params = &Params;

	DBS.LoadDB(QFN);
	DBS.Setup();
	StatSig::Init(DBS.GetDBChainCount());

	OpenOutputFiles();
	ChainReader2 CR;
	CR.Open(DBFN);
	DBS.RunQuery(CR);
	CloseOutputFiles();
	}

void cmd_search()
	{
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

	DSSParams Params;
	Params.SetDSSParams(DM_UseCommandLineOption);
	const string &PatternStr = prefiltermu_pattern;
	asserta(PatternStr == "1110011");

	string MuFilterTsvFN;
	GetTmpFileName(MuFilterTsvFN);
	Log("MuFilterTsvFN=%s\n", MuFilterTsvFN.c_str());

	MuSeqSource QSS;
	MuSeqSource DBSS;

	QSS.OpenChains(QueryFN, Params);

	if (optset_dbmu)
		DBSS.OpenFasta(opt(dbmu));
	else
		DBSS.OpenChains(DBFN, Params);

	SeqDB MuQueryDB;
	MuQueryDB.FromSS(QSS);

	uint DBSize = MuPreFilter(Params, MuQueryDB, DBSS, MuFilterTsvFN);
	StatSig::Init(DBSize);

	DSSParams Params2;
	Params2.SetDSSParams(DM_AlwaysSensitive);
	PostMuFilter(Params2, MuFilterTsvFN, QueryFN, DBFN, opt(output));

	if (!opt(keeptmp))
		DeleteStdioFile(MuFilterTsvFN);
	}
