#include "myutils.h"

/***
950K Jul  6 08:17 rdrp_pdb70		# hits in TSV format with a few null bytes inserted (not understood)
 597 Jul  6 08:17 rdrp_pdb70.index	# TSV file giving start byte position and record length in bytes for each hit
   4 Jul  6 08:17 rdrp_pdb70.dbtype	# 4-byte file always contains hex 00C000

# head *.index
0       0       5472
1       5472    48042
2       53514   12741
3       66255   52115

# tail *.index
35      897583  47799
36      945382  3833
37      949215  19512
38      968727  3284
***/

void cmd_mmseqs_index_dump()
	{
	const string &Prefix = g_Arg1;
	const string &HitsFN = Prefix;
	const string &IndexFN = Prefix + ".index";
	const string &DbtypeFN = Prefix + ".dbtype";

	FILE *fout = 0;
	if (optset_output)
		fout = CreateStdioFile(opt(output));

	FILE *f = OpenStdioFile(DbtypeFN);
	uint64_t Dbtype_filesize = GetStdioFileSize64(f);
	asserta(Dbtype_filesize == 4);
	uint32_t u;
	ReadStdioFile(f, &u, 4);
	CloseStdioFile(f);
	ProgressLog("0x%04x  %s\n", u, DbtypeFN.c_str());

	FILE *fhits = OpenStdioFile(HitsFN);
	FILE *findex = OpenStdioFile(IndexFN);
	uint64_t hits_filesize = GetStdioFileSize64(f);
	string line;
	vector<string> flds;
	uint recnr = 0;
	uint nextrecpos = 0;
	uint nonprint = 0;
	uint hitcount = 0;
	while (ReadLineStdioFile(findex, line))
		{
		Split(line, flds, '\t');
		asserta(flds.size() == 3);
		uint recidx = StrToUint(flds[0]);
		asserta(recidx == recnr);
		++recnr;
		uint recpos = StrToUint(flds[1]);
		uint reclen = StrToUint(flds[2]);
		asserta(recpos == nextrecpos);
		nextrecpos += reclen;
		asserta(reclen > 0);
		byte *buffer = myalloc(byte, reclen);
		ReadStdioFile(fhits, recpos, buffer, reclen);
		asserta(buffer[reclen-1] == 0);
		if (fout != 0)
			{
			fprintf(fout, "index\t%u\t%u\n", recpos, reclen);
			for (uint i = 0; i < reclen-1; ++i)
				{
				char c = buffer[i];
				if (c == '\n')
					{
					fprintf(fout, "\n");
					++hitcount;
					}
				else
					{
					if (!isprint(c) && c != '\t')
						{
						++nonprint;
						fprintf(fout, "@");
						}
					else
						fprintf(fout, "%c", c);
					}
				}
			fprintf(fout, "\n");
			}
		}
	asserta(nextrecpos == hits_filesize);

	CloseStdioFile(fhits);
	CloseStdioFile(findex);
	CloseStdioFile(fout);
	ProgressLog("%u records, %u hits, %u non-printing bytes\n",
				recnr, hitcount, nonprint);
	}
