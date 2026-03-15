#include "myutils.h"

void trunc_label(string &Label);

void cmd_make_pairs_binary()
	{
	asserta(optset_output);
	const string &lookupfn = g_Arg1;
	const string &hitsfn = opt(input);	// first two fields are query,target

	vector<string> lines;
	ReadLinesFromFile(lookupfn, lines);

	vector<string> flds;
	map<string, uint> label2idx;
	for (uint idx = 0; idx < SIZE(lines); ++idx)
		{
		Split(lines[idx], flds, '\t');
		asserta(SIZE(flds) == 2);
		string label = flds[0];
		trunc_label(label);
		label2idx[label] = idx;
		}

	FILE *f = OpenStdioFile(hitsfn);
	FILE *fOut = CreateStdioFile(opt(output));
	string line;
	uint npair = 0;
	uint missing = 0;
	uint bytes = 0;
	while (ReadLineStdioFile(f, line))
		{
		Split(line, flds, '\t');
		asserta(SIZE(flds) >= 2);
		string q = flds[0];
		string t = flds[1];
		trunc_label(q);
		trunc_label(t);
		map<string, uint>::const_iterator iterq = label2idx.find(q);
		map<string, uint>::const_iterator itert = label2idx.find(t);
		if (iterq == label2idx.end() || itert == label2idx.end())
			{
			++missing;
			continue;
			}
		++npair;
		uint16_t idxq = uint16_t(iterq->second);
		uint16_t idxt = uint16_t(iterq->second);
		WriteStdioFile(fOut, &idxq, sizeof(idxq));
		WriteStdioFile(fOut, &idxt, sizeof(idxt));
		bytes += 2*sizeof(uint16_t);
		}
	CloseStdioFile(f);
	CloseStdioFile(fOut);
	ProgressLog("%u pairs, %u missing, %s bytes\n",
		npair, missing, MemBytesToStr(double(bytes)));
	}