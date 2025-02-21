#include "myutils.h"
#include "dss.h"
#include "seqdb.h"
#include "museqsource.h"
#include "dbsearcher.h"
#include "search.h"
#include "output.h"

uint MuPreFilter(const DSSParams &Params,
			  SeqDB &QueryDB,
			  MuSeqSource &FSS,
			  const string &OutputFN);

void SelfSearch()
	{
	const string &QFN = g_Arg1;
	if (optset_db)
		Die("-db not used for -selfsearch");

	DBSearcher DBS;
	DSSParams Params;
	Params.SetDSSParams(DM_UseCommandLineOption, SCOP40_DBSIZE);
	DBS.m_Params = &Params;

	DBS.LoadDB(QFN);
	DBS.Setup();
	Params.m_DBSize = (float) DBS.GetDBSize();
	if (optset_dbsize)
		Params.m_DBSize = (float) opt_dbsize;

	OpenOutputFiles();
	DBS.RunSelf();
	CloseOutputFiles();
	}

static void Search_NoMuFilter()
	{
	if (!optset_db)
		Die("-db required");

	const string &QFN = g_Arg1;
	const string &DBFN = opt_db;

	DBSearcher DBS;
	DSSParams Params;
	Params.SetDSSParams(DM_UseCommandLineOption, SCOP40_DBSIZE);
	DBS.m_Params = &Params;

	DBS.LoadDB(QFN);
	DBS.Setup();

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
	const string &DBFN = string(opt_db);

	if (!EndsWith(DBFN, ".bca"))
		Die(".bca format required for -db");

	DSSParams Params;
	Params.SetDSSParams(DM_UseCommandLineOption, SCOP40_DBSIZE);
	const string &PatternStr = Params.m_MuPrefPatternStr;
	asserta(PatternStr == "1110011");

	string MuFilterTsvFN;
	GetTmpFileName(MuFilterTsvFN);
	Log("MuFilterTsvFN=%s\n", MuFilterTsvFN.c_str());

	MuSeqSource QSS;
	MuSeqSource DBSS;
	QSS.Open(QueryFN, Params);
	DBSS.Open(DBFN, Params);

	SeqDB MuQueryDB;
	MuQueryDB.FromSS(QSS);

	float MaxEvalue = 10;
	if (optset_evalue)
		MaxEvalue = (float) opt_evalue;

	uint DBSize = MuPreFilter(Params, MuQueryDB, DBSS, MuFilterTsvFN);

	DSSParams Params2;
	Params2.SetDSSParams(DM_AlwaysSensitive, SCOP40_DBSIZE);
	PostMuFilter(Params2, MuFilterTsvFN, QueryFN, DBFN, MaxEvalue, opt_output);

	if (!opt_keeptmp)
		DeleteStdioFile(MuFilterTsvFN);
	}
