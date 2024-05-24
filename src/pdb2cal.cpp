#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"

void cmd_pdb2cal()
	{
	ChainReader CR;
	CR.Open(g_Arg1);

	FILE *fOut = CreateStdioFile(opt_output);

	PDBChain Chain;
	uint Count = 0;
	for (;;)
		{
		bool Ok = CR.GetNext(Chain);
		if (!Ok)
			break;

		if (++Count%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u converted >%s\r",
			  sPct.c_str(), Count, Chain.m_Label.c_str());
			}
		Chain.ToCal(fOut);
		}
	Progress("100.0%% done, %u chains converted\n", Count);
	CloseStdioFile(fOut);
	}
