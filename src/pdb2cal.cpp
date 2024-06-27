#include "myutils.h"
#include "pdbchain.h"
#include "chainreader2.h"

void cmd_pdb2cal()
	{
	PDBFileScanner FS;
	FS.Open(g_Arg1);

	ChainReader2 CR;
	CR.Open(FS);

	FILE *fOut = CreateStdioFile(opt_output);

	uint Count = 0;
	for (;;)
		{
		PDBChain *ptrChain = CR.GetNext();
		if (ptrChain == 0)
			break;
		PDBChain &Chain = *ptrChain;
		Chain.ToCal(fOut);
		delete ptrChain;
		}
	Progress("100.0%% done, %u chains converted\n", Count);
	CloseStdioFile(fOut);
	}

void cmd_pdb2fasta()
	{
	PDBFileScanner FS;
	FS.Open(g_Arg1);

	ChainReader2 CR;
	CR.Open(FS);

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
