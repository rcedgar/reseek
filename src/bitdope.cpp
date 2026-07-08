#include "myutils.h"
#include "triangle.h"
#include "lookup.h"
#include "bitdope.h"

void trunc_label(string &Label);

uint8_t *read_bitdope(const string &fn,
	uint32_t &ndom, uint32_t &nhit)
	{
	bitdope dope;
	dope.from_file(fn);
	ndom = dope.m_ndom;
	nhit = dope.m_nhit;
	return dope.m_dope;
	}

#if 0
void cmd_bitdope_stats()
	{
	bitdope dope;
	dope.from_file(g_Arg1);
	ProgressLog("ndom=%u  nhit=%u  %s\n",
		dope.m_ndom, dope.m_nhit, g_Arg1.c_str());
	}
#endif

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
	uint ntp = 0;
	uint discarded = 0;
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
		if (look.is_tp_ij(minidx, maxidx))
			++ntp;
		else
			{
			if (optset_truth)
				{
				++discarded;
				continue;
				}
			}
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
	ProgressLog("nhit %u, ntp %u\n", nhit, ntp);
	if (nbit != nhit)
		Die("nbit %u, nhit %u", nbit, nhit);
	if (discarded > 0)
		ProgressLog("discarded %u\n", discarded);

	FILE *fOut = CreateStdioFile(opt(output));
	WriteStdioFile(fOut, &MAGIC, sizeof(MAGIC));
	WriteStdioFile(fOut, &ndom, sizeof(ndom));
	WriteStdioFile(fOut, bitvec, bytes);
	WriteStdioFile(fOut, &MAGIC, sizeof(MAGIC));
	CloseStdioFile(fOut);

	bitdope dope2;
	dope2.m_look = &look;
	dope2.from_file(opt(output));
	asserta(dope2.m_ndom == ndom);
	if (dope2.m_nhit != nhit)
		Die("nhit2 %u, nhit %u", dope2.m_nhit, nhit);
	}
