#include "myutils.h"
#include "dbsearcher.h"
#include "mx.h"
#include "pdbchain.h"
#include "xdpmem.h"
#include "cigar.h"
#include "timing.h"
#include <thread>

void cmd_search()
	{
	const string &QFN = g_Arg1;
	const string &DBFN = opt_db;

	DBSearcher DBS;
	DSSParams Params;
	Params.SetFromCmdLine(10000);
	DBS.m_Params = &Params;

	bool Self = (DBFN == "");

	if (Self)
		DBS.LoadDB(DBFN);
	else
		DBS.LoadDB(QFN);

	Params.m_DBSize = (float) DBS.GetDBSize();
	DBS.Setup();

	DBS.m_fTsv = CreateStdioFile(opt_output);
	DBS.m_fAln = CreateStdioFile(opt_aln);
	DBS.m_fFasta2 = CreateStdioFile(opt_fasta2);
	ResetTimers();
	if (Self)
		DBS.RunSelf();
	else
		DBS.RunSearch();
	CloseStdioFile(DBS.m_fTsv);
	}
