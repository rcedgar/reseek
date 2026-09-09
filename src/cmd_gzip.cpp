#include "myutils.h"

void cmd_gunzip()
	{
	const string &InputFileName = g_Arg1;
	const string &OutputFileName = opt(output);

	FILE *fIn = OpenGzipFile(InputFileName);
	FILE *fOut = CreateStdioFile(OutputFileName);

	const unsigned M = 1024*1024;
	byte *Buffer = myalloc(byte, M);

	for (;;)
		{
		uint32 n = ReadGzipFile(fIn, Buffer, M);
		if (n == 0)
			break;
		WriteStdioFile(fOut, Buffer, n);
		}

	CloseGzipFile(fIn);
	CloseStdioFile(fOut);
	}

void cmd_gunzip_lines()
	{
	const string &InputFileName = g_Arg1;

	vector<string> Lines;
	ReadLinesFromGzipFile(InputFileName, Lines);

	if (!optset_output)
		return;

	FILE *fOut = CreateStdioFile(opt(output));
	for (uint i = 0; i < SIZE(Lines); ++i)
		{
		fputs(Lines[i].c_str(), fOut);
		fputc('\n', fOut);
		}
	CloseStdioFile(fOut);
	}
