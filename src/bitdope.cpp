#include "myutils.h"
#include "triangle.h"
#include "lookup.h"

void trunc_label(string &Label);

static const uint32_t MAGIC	= 0xd05e;

/***
Read TSV 1=query 2=target from mufilter.
Save sorted list of 7-character used labels (SCOP domains)
Make bitvector with pairs which passed.

0123456
d2eyqa6
d7reqb1
***/

uint8_t *read_bitdope(const string &fn, uint32_t &ndom, uint32_t &nhit)
	{
	uint32_t magic;
	FILE *f = OpenStdioFile(fn);
	ReadStdioFile(f, &magic, sizeof(magic));
	asserta(magic == MAGIC);
	ReadStdioFile(f, &ndom, sizeof(ndom));
	uint32_t K = triangle_get_K(ndom);
	uint32_t bytes = (K + 7)/8;
	uint8_t *bitvec = myalloc(uint8_t, bytes);

	ReadStdioFile(f, bitvec, bytes);
	ReadStdioFile(f, &magic, sizeof(magic));
	asserta(magic == MAGIC);
	CloseStdioFile(f);

	for (uint i = 0; i < bytes; ++i)
		{
		uint8_t b = bitvec[i];
		for (uint j = 0; j < 8; ++j)
			{
			if (b & (1 << j))
				++nhit;
			}
		}
	return bitvec;
	}

void cmd_bitdope()
	{
	asserta(optset_output);
	asserta(optset_lookup);
	const string &hitsfn = g_Arg1;

	lookup look;
	look.from_tsv(opt(lookup));
	const uint32_t ndom = look.get_ndom();
	uint32_t K = triangle_get_K(ndom);

	uint32_t bytes = (K + 7)/8;
	uint8_t *bitvec = myalloc(uint8_t, bytes);
	memset(bitvec, 0, bytes);

	FILE *f = OpenStdioFile(hitsfn);
	string line;
	vector<string> flds;
	uint nhit = 0;
	while (ReadLineStdioFile(f, line))
		{
		Split(line, flds, '\t');
		asserta(flds.size() == 2);
		const string &labelq = flds[0];
		const string &labelt = flds[1];
		uint idxq = look.get_domidx(labelq, true);
		uint idxt = look.get_domidx(labelt, true);
		if (idxq == UINT_MAX || idxt == UINT_MAX)
			continue;
		uint minidx = min(idxq, idxt);
		uint maxidx = max(idxq, idxt);
		uint k = triangle_ij_to_k(minidx, maxidx, ndom);
		bitvec[k/8] |= (1 << (k%8));
		++nhit;
		}
	CloseStdioFile(f);

	uint nbit = 0;
	for (uint i = 0; i < bytes; ++i)
		{
		uint8_t b = bitvec[i];
		for (uint j = 0; j < 8; ++j)
			{
			if (b & (1 << j))
				++nbit;
			}
		}
	ProgressLog("nbit %u\n", nbit);
	asserta(nbit == nhit);

	FILE *fOut = CreateStdioFile(opt(output));
	WriteStdioFile(fOut, &MAGIC, sizeof(MAGIC));
	WriteStdioFile(fOut, &ndom, sizeof(ndom));
	WriteStdioFile(fOut, bitvec, bytes);
	WriteStdioFile(fOut, &MAGIC, sizeof(MAGIC));
	CloseStdioFile(fOut);

	uint32_t ndom2, nhit2;
	read_bitdope(opt(output), ndom2, nhit2);
	asserta(ndom2 == ndom);
	asserta(nhit2 == nhit);
	}
