#include "myutils.h"
#include "dbsearcher.h"
#include "mx.h"
#include "pdbchain.h"
#include "xdpmem.h"
// #include "outputfiles.h"
#include "cigar.h"
#include "timing.h"
#include <thread>

void cmd_search()
	{
	const string &QCalFN = g_Arg1;
	const string &DBFN = opt_db;
	DBSearcher DBS;
	DSSParams Params;
	Params.SetFromCmdLine();
	DBS.ReadChains(QCalFN, DBFN);
	Params.m_DBSize = (float) DBS.GetDBSize();
	if (optset_dbsize)
		Params.m_DBSize = (float) opt_dbsize;

	DBS.Setup(Params);
	DBS.m_MaxEvalue = FLT_MAX;
	DBS.m_fTsv = CreateStdioFile(opt_output);
	DBS.m_fAln = CreateStdioFile(opt_aln);
	ResetTimers();
	DBS.Run();
	DSSAligner::Stats();
	CloseStdioFile(DBS.m_fTsv);
	}
