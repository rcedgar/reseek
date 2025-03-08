#include "myutils.h"
#include "museqsource.h"
#include "pdbchain.h"
#include "dss.h"
#include "seqinfo.h"

void cmd_convert2mu()
	{
	optset_fast = true;
	opt(fast) = true;

	uint MinChainLength = 1;
	if (optset_minchainlength)
		MinChainLength = opt(minchainlength);

	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast, SCOP40_DBSIZE);
	MuSeqSource MSS;
	MSS.OpenChains(g_Arg1, Params);

	if (optset_output)
		{
		opt(fasta) = opt(output);
		optset_fasta = true;
		}

	FILE *fFasta = CreateStdioFile(opt(fasta));

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
