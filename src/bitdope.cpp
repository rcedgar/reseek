#include "myutils.h"
#include "triangle.h"

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

uint8_t *read_bitdope(
	const string &fn,
	vector<string> &labels)
	{
	uint32_t magic, nlab;
	FILE *f = OpenStdioFile(fn);
	ReadStdioFile(f, &magic, sizeof(magic));
	asserta(magic == MAGIC);
	ReadStdioFile(f, &nlab, sizeof(nlab));
	uint32_t K = triangle_get_K(nlab);
	uint32_t bytes = (K + 7)/8;
	uint8_t *bitvec = myalloc(uint8_t, bytes);

	char label[8];
	for (uint i = 0; i < nlab; ++i)
		{
		ReadStdioFile(f, label, 8);
		string slabel = string(label);
		labels.push_back(label);
		}
	ReadStdioFile(f, bitvec, bytes);
	ReadStdioFile(f, &magic, sizeof(magic));
	asserta(magic == MAGIC);
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
	ProgressLog("read_bitdope() nbit %u\n", nbit);
	return bitvec;
	}

static void add_label(
	const string &label,
	unordered_map<string, uint> &label2idx,
	vector<string> &labels)
	{
	unordered_map<string, uint>::const_iterator iterq =
			label2idx.find(label);
	if (iterq == label2idx.end())
		{
		uint idx = uint(labels.size());
		label2idx[label] = idx;
		labels.push_back(label);
		}
	}

void cmd_bitdope()
	{
	asserta(optset_output);
	const string &fn = g_Arg1;

	FILE *f = OpenStdioFile(fn);
	string line;
	vector<string> flds;
	unordered_map<string, uint> label2idx;
	vector<string> labelqs;
	vector<string> labelts;
	vector<string> uniq_labels;
	while (ReadLineStdioFile(f, line))
		{
		Split(line, flds, '\t');
		asserta(flds.size() == 2);
		string labelq = flds[0];
		string labelt = flds[1];
		trunc_label(labelq);
		trunc_label(labelt);
		if (labelq == labelt)
			continue;
		asserta(strlen(labelq.c_str()) == 7);
		asserta(strlen(labelt.c_str()) == 7);

		labelqs.push_back(labelq);
		labelts.push_back(labelt);

		add_label(labelq, label2idx, uniq_labels);
		add_label(labelt, label2idx, uniq_labels);
		}
	CloseStdioFile(f);

	const uint32_t nlab = uint(uniq_labels.size());
	uint32_t K = triangle_get_K(nlab);

	uint32_t bytes = (K + 7)/8;
	uint8_t *bitvec = myalloc(uint8_t, bytes);
	memset(bitvec, 0, bytes);

	const size_t nhit = labelqs.size();
	asserta(labelts.size() == nhit);

	ProgressLog("%u labels, K %u, hits %u (%.1f%%)\n",
		nlab, K, nhit, GetPct(double(nhit), double(K)));
	for (size_t i = 0; i < nhit; ++i)
		{
		const string &labelq = labelqs[i];
		const string &labelt = labelts[i];
		uint idxq = label2idx[labelq];
		uint idxt = label2idx[labelt];
		uint minidx = min(idxq, idxt);
		uint maxidx = max(idxq, idxt);
		uint k = triangle_ij_to_k(minidx, maxidx, nlab);

		bitvec[k/8] |= (1 << (k%8));
		}

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

	FILE *fOut = CreateStdioFile(opt(output));

	WriteStdioFile(fOut, &MAGIC, sizeof(MAGIC));
	WriteStdioFile(fOut, &nlab, sizeof(nlab));
	for (uint i = 0; i < nlab; ++i)
		{
		const string &label = uniq_labels[i];
		asserta(strlen(label.c_str()) == 7);
		WriteStdioFile(fOut, label.c_str(), 8);
		}
	WriteStdioFile(fOut, bitvec, bytes);
	WriteStdioFile(fOut, &MAGIC, sizeof(MAGIC));
	CloseStdioFile(fOut);

	vector<string> labels2;
	read_bitdope(opt(output), labels2);
	asserta(labels2 == uniq_labels);
	ProgressLog("label check passed\n");
	}
