#include "myutils.h"
#include "dbsearcher.h"
#include "chainreader2.h"
#include "output.h"

void cmd_smallsearch()
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
