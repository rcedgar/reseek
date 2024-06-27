#include "myutils.h"
#include "chainreader2.h"

void ReadChains(const string &FileName, vector<PDBChain *> &Chains)
	{
	PDBFileScanner FS;
	FS.Open(FileName);

	ChainReader2 CR;
	CR.Open(FS);
	for (;;)
		{
		PDBChain *Chain = CR.GetNext();
		if (Chain == 0)
			break;
		Chains.push_back(Chain);
		}
	}
