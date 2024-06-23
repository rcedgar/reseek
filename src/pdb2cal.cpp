#include "myutils.h"
#include "pdbchain.h"
#include "chainreader2.h"

void cmd_pdb2cal()
	{
	ChainReader2 CR;
	CR.Open(g_Arg1);

	FILE *fOut = CreateStdioFile(opt_output);

	uint Count = 0;
	for (;;)
		{
		PDBChain *ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		PDBChain &Chain = *ptrChain;

		if (++Count%100 == 0)
			{
			Progress("%u converted >%s\r",
			  Count, Chain.m_Label.c_str());
			}
		Chain.ToCal(fOut);
		}
	Progress("100.0%% done, %u chains converted\n", Count);
	CloseStdioFile(fOut);
	}

void cmd_pdb2fasta()
	{
	ChainReader2 CR;
	CR.Open(g_Arg1);

	FILE *fOut = CreateStdioFile(opt_output);

	uint Count = 0;
	for (;;)
		{
		PDBChain *ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		PDBChain &Chain = *ptrChain;
		SeqToFasta(fOut, Chain.m_Label, Chain.m_Seq);
		delete ptrChain;
		}
	Progress("%u chains converted\n", Count);
	CloseStdioFile(fOut);
	}
