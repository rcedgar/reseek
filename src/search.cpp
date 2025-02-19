#include "myutils.h"
#include "dss.h"
#include "seqdb.h"
#include "museqsource.h"
#include "dbsearcher.h"
#include "search.h"
#include "output.h"

void SelfSearch()
	{
	const string &QFN = g_Arg1;
	if (optset_db)
		Die("-db not used for -selfsearch");

	DBSearcher DBS;
	DSSParams Params;
	Params.SetFromCmdLine(10000);
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
	Params.SetFromCmdLine(10000);
	DBS.m_Params = &Params;

	DBS.LoadDB(QFN);

	if (optset_dbsize)
		Params.m_DBSize = (float) opt_dbsize;
	else
		{
		Warning("-dbsize not set, defaulting to effective size 10,000 chains");
		Params.m_DBSize = 10000;
		}
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
	Params.SetFromCmdLine(10000);
	const string &PatternStr = Params.m_PatternStr;
	asserta(PatternStr == "111");

	string MuFilterTsvFN;
	GetTmpFileName(MuFilterTsvFN);
	Log("MuFilterTsvFN=%s\n", MuFilterTsvFN.c_str());

	MuSeqSource QSS;
	MuSeqSource DBSS;
	QSS.Open(QueryFN, Params);
	if (optset_input2)
		{
		Progress("open %s\n", opt_input2);
		DBSS.OpenFasta(opt_input2);
		}
	else
		DBSS.Open(DBFN, Params);

	SeqDB MuQueryDB;
	MuQueryDB.FromSS(QSS);

	float MaxEvalue = 10;
	if (optset_evalue)
		MaxEvalue = (float) opt_evalue;

	Die("TODO -- MuFilter");
#if 0
	uint DBSize = MuFilter(Params, MuQueryDB, DBSS, MuFilterTsvFN);
	if (optset_dbsize)
		DBSize = uint(opt_dbsize);
	Params.m_DBSize = float(DBSize);
#endif
	PostMuFilter(Params, MuFilterTsvFN, QueryFN, DBFN, MaxEvalue, opt_output);

	if (!opt_keeptmp)
		DeleteStdioFile(MuFilterTsvFN);
	}
