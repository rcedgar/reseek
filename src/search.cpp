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
	const string &QCalFN = g_Arg1;
	const string &DBFN = opt_db;
	DBSearcher DBS;
	DSSParams Params;
	DBS.ReadChains(QCalFN, DBFN);
	uint DBSize = DBS.GetDBSize();
	Params.SetFromCmdLine(DBS.GetDBSize());

	DBS.Setup(Params);
	DBS.m_fTsv = CreateStdioFile(opt_output);
	DBS.m_fAln = CreateStdioFile(opt_aln);
	DBS.m_fFasta2 = CreateStdioFile(opt_fasta2);
	ResetTimers();
	DBS.Run();
	CloseStdioFile(DBS.m_fTsv);
	}
