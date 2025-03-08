#include "myutils.h"
#include "output.h"

FILE *g_fTsv;
FILE *g_fAln;
FILE *g_fFasta2;

void OpenOutputFiles()
	{
	g_fTsv = CreateStdioFile(opt(output));
	g_fAln = CreateStdioFile(opt(aln));
	g_fFasta2 = CreateStdioFile(opt(fasta2));
	}

void CloseOutputFiles()
	{
	CloseStdioFile(g_fTsv);
	CloseStdioFile(g_fAln);
	CloseStdioFile(g_fFasta2);
	}
