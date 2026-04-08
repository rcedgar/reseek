#include "myutils.h"
#include "triangle.h"
#include "lookup.h"

void trunc_label(string &Label);

static const uint32_t MAGIC	= 0xd05e;

uint8_t *read_bitdope(const string &fn,
	uint32_t &ndom, uint32_t &nhit)
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

	nhit = 0;
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

/***
The -bitdope command creates a bit-vector file representing an 
all-vs-all search. A bit is 1/0 if hit is/not present.
Only the upper triangle is represented, i.e. only pairs (i,j) where
	0 <= i <= j < N.
The diagonal (i,i) is included but is never used because these
represent trivial self-hits (basically I was lazy and didn't
update triangle.h to exclude the diagonal).
***/
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
		asserta(flds.size() >= 2);
		const string &labelq = flds[0];
		const string &labelt = flds[1];
		uint idxq = look.get_domidx(labelq, true);
		uint idxt = look.get_domidx(labelt, true);
		if (idxq == UINT_MAX || idxt == UINT_MAX)
			continue;
		uint minidx = min(idxq, idxt);
		uint maxidx = max(idxq, idxt);
		uint k = triangle_ij_to_k(minidx, maxidx, ndom);
		const uint8_t thebit = (1 << (k%8));
		if ((bitvec[k/8] & thebit) == 0)
			{
			bitvec[k/8] |= thebit;
			++nhit;
			}
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
	ProgressLog("nhit %u\n", nhit);
	if (nbit != nhit)
		Die("nbit %u, nhit %u", nbit, nhit);

	FILE *fOut = CreateStdioFile(opt(output));
	WriteStdioFile(fOut, &MAGIC, sizeof(MAGIC));
	WriteStdioFile(fOut, &ndom, sizeof(ndom));
	WriteStdioFile(fOut, bitvec, bytes);
	WriteStdioFile(fOut, &MAGIC, sizeof(MAGIC));
	CloseStdioFile(fOut);

	uint32_t ndom2, nhit2;
	const uint8_t *bitvec2 =
		read_bitdope(opt(output), ndom2, nhit2);
	asserta(ndom2 == ndom);
	if (nhit2 != nhit)
		Die("nhit2 %u, nhit %u", nhit2, nhit);
	}
