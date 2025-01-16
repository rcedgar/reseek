#include "myutils.h"
#include "dss.h"
#include "seqdb.h"
#include "museqsource.h"
#include "search.h"

void cmd_search()
	{
	if (!optset_db)
		Die("-db option required");

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

	uint DBSize = MuFilter(Params, MuQueryDB, DBSS, MuFilterTsvFN);
	if (optset_dbsize)
		DBSize = uint(opt_dbsize);
	Params.m_DBSize = float(DBSize);

	PostMuFilter(Params, MuFilterTsvFN, QueryFN, DBFN, MaxEvalue, opt_output);

	if (!opt_keeptmp)
		DeleteStdioFile(MuFilterTsvFN);
	}
