#include "myutils.h"
#include "dbsearcher.h"
#include "output.h"

void cmd_selfsearch()
	{
	const string &QFN = g_Arg1;
	if (optset_db)
		Die("-db not used for -selfsearch");

	DBSearcher DBS;
	DSSParams Params;
	Params.SetFromCmdLine(10000);
	DBS.m_Params = &Params;

	DBS.LoadDB(QFN);
	Params.m_DBSize = (float) DBS.GetDBSize();
	if (optset_dbsize)
		Params.m_DBSize = (float) opt_dbsize;

	OpenOutputFiles();
	DBS.RunSelf();
	CloseOutputFiles();
	}
