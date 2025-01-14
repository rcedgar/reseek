#include "myutils.h"
#include "museqsource.h"
#include "pdbchain.h"
#include "dss.h"
#include "seqinfo.h"

void cmd_convert2mu()
	{
	optset_fast = true;
	opt_fast = true;

	uint MinChainLength = 1;
	if (optset_minchainlength)
		MinChainLength = opt_minchainlength;

	DSSParams Params;
	Params.SetFromCmdLine(10000);
	MuSeqSource MSS;
	MSS.Open(g_Arg1, Params);

	FILE *fFasta = CreateStdioFile(opt_fasta);

	ObjMgr OM;
	uint n = 0;
	time_t t_prog = time(0);
	SeqInfo *SI = OM.GetSeqInfo();
	for (;;)
		{
		if (n%100 == 0)
			{
			time_t now = time(0);
			if (now - t_prog> 0)
				{
				t_prog = now;
				Progress("%u chains converted  (%s)   \r", n, IntToStr(n));
				}
			}
		bool Ok = MSS.GetNext(SI);
		if (!Ok)
			break;
		SeqToFasta(fFasta, SI->m_Label, (const char *) SI->m_Seq, SI->m_L);
		++n;
		}

	Progress("%u chains converted  (%s)   \n", n, IntToStr(n));
	CloseStdioFile(fFasta);
	}
