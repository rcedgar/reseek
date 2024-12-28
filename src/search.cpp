#include "myutils.h"
#include "dbsearcher.h"
#include "chainreader2.h"
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
	DBS.LoadDB(QFN);

	if (Self)
		Params.m_DBSize = (float) DBS.GetDBSize();
	else
		{
		if (optset_dbsize)
			Params.m_DBSize = (float) opt_dbsize;
		else
			{
			Warning("-dbsize not set, defaulting to effective size 10,000 chains");
			Params.m_DBSize = 10000;
			}
		}
	DBS.Setup();

	DBS.m_fTsv = CreateStdioFile(opt_output);
	DBS.m_fAln = CreateStdioFile(opt_aln);
	DBS.m_fFasta2 = CreateStdioFile(opt_fasta2);
	ResetTimers();
	if (Self)
		DBS.RunSelf();
	else
		{
		ChainReader2 CR;
		CR.Open(DBFN);
		DBS.RunQuery(CR);
		}
	CloseStdioFile(DBS.m_fTsv);
	}
