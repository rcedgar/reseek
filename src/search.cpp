#include "myutils.h"
#include "dss.h"
#include "seqdb.h"
#include "museqsource.h"
#include "dbsearcher.h"
#include "output.h"
#include "statsig.h"
#include "prefiltermuparams.h"

uint MuPreFilter(SeqDB &QueryDB,
			  MuSeqSource &FSS,
			  const string &OutputFN);

void PostMuFilter(const string &MuFilterTsvFN,
				  const string &QueryCAFN,
				  const string &DBBCAFN,
				  const string &HitsFN);

void SelfSearch()
	{
	const string &QFN = g_Arg1;
	if (optset_db)
		Die("-db not used for -selfsearch");

	DBSearcher DBS;

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

	string MuFilterTsvFN;
	GetTmpFileName(MuFilterTsvFN);
	Log("MuFilterTsvFN=%s\n", MuFilterTsvFN.c_str());

	MuSeqSource QSS;
	MuSeqSource DBSS;

	QSS.OpenChains(QueryFN);

	if (optset_dbmu)
		DBSS.OpenFasta(opt(dbmu));
	else
		DBSS.OpenChains(DBFN);

	SeqDB MuQueryDB;
	MuQueryDB.FromSS(QSS);

	uint DBSize = MuPreFilter(MuQueryDB, DBSS, MuFilterTsvFN);
	StatSig::Init(DBSize);

	DSSParams::ReInit(DM_AlwaysSensitive);
	PostMuFilter(MuFilterTsvFN, QueryFN, DBFN, opt(output));

	if (!opt(keeptmp))
		DeleteStdioFile(MuFilterTsvFN);
	}
