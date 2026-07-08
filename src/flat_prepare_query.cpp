#include "myutils.h"
#include "flat_chain.h"
#include "flat_params.h"
#include "bcadata.h"
#include "chaq.h"
#include "xdpmem.h"
#include <set>

float ViterbiFastMem_Blosum62(XDPMem &Mem, const char *A, uint LA,
	const char *B, uint LB, string &Path);

static double GetPctId(const string &Seq_i, const string &Seq_j)
	{
	const uint LA = SIZE(Seq_i);
	const uint LB = SIZE(Seq_j);
	if (LA == LB && Seq_i == Seq_j)
		return 100;
	const char *A = Seq_i.c_str();
	const char *B = Seq_j.c_str();
	string Path;
	XDPMem Mem;
	ViterbiFastMem_Blosum62(Mem, A, LA, B, LB, Path);
	const uint ColCount = SIZE(Path);
	uint PosA = 0;
	uint PosB = 0;
	uint Ids = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			if (Seq_i[PosA] == Seq_j[PosB])
				++Ids;
			++PosA;
			++PosB;
			break;
		case 'D':
			++PosA;
			break;
		case 'I':
			++PosB;
			break;
			}
		}
	double PctId = (100.0*Ids)/ColCount;
	return PctId;
	}

static string GetFlatChainSeq(const flat_chain_t &chain)
	{
	const uint L = chain.get_length();
	string seq;
	seq.reserve(L);
	for (uint pos = 0; pos < L; ++pos)
		seq.push_back(chain.get_aa(pos));
	return seq;
	}

void cmd_flat_prepare_query()
	{
	if (!optset_bcb)
		Die("Must specify -bcb");

	BCAData BCA;
	BCA.Create(opt(bcb), true);

	FILE *fOut = CreateStdioFile(opt(output));

	vector<flat_chain_t *> input_chains;
	read_flat_chains(g_Arg1, input_chains);
	const uint input_chain_count = uint(input_chains.size());

	const double min_pct_id = 90;
	const uint min_len = (optset_minchainlength ? opt(minchainlength) : 1);
	const uint max_chains = (optset_n ? opt(n) : 4);

	vector<flat_chain_t *> output_chains;
	set<uint> deleted_chain_idxs;
	uint too_short = 0;
	uint nr_queries = 0;

	for (uint i = 0; i < input_chain_count; ++i)
		{
		if (deleted_chain_idxs.find(i) != deleted_chain_idxs.end())
			continue;

		const flat_chain_t &chain_i = *input_chains[i];
		const string &label_i = chain_i.m_label;
		const string seq_i = GetFlatChainSeq(chain_i);
		const uint L = uint(seq_i.size());

		Pf(fOut, "%u\t%s\t%u", i, label_i.c_str(), L);
		if (L < min_len)
			{
			Pf(fOut, "\tshort\n");
			++too_short;
			continue;
			}

		if (L > flat_params::m_maxL)
			{
			Pf(fOut, "\ttoolong\n");
			continue;
			}

		bool del = false;
		if (nr_queries >= max_chains)
			{
			del = true;
			Pf(fOut, "\ttoomany\n");
			continue;
			}

		for (uint j = 0; j < i; ++j)
			{
			if (deleted_chain_idxs.find(j) != deleted_chain_idxs.end())
				continue;
			const flat_chain_t &chain_j = *input_chains[j];
			const string seq_j = GetFlatChainSeq(chain_j);
			if (SIZE(seq_j) < min_len)
				continue;
			double pct_id = GetPctId(seq_i, seq_j);
			if (pct_id >= min_pct_id)
				{
				Pf(fOut, "\t%.1f%%%u\n", pct_id, j);
				deleted_chain_idxs.insert(i);
				del = true;
				break;
				}
			}

		if (!del)
			{
			output_chains.push_back(input_chains[i]);
			++nr_queries;
			Pf(fOut, "\tquery\n");
			}
		}

	chaq_vecs2 cv;
	chaq::alloc_chaq_vecs2(cv, flat_params::m_maxL);
	const uint output_chain_count = uint(output_chains.size());
	for (uint i = 0; i < output_chain_count; ++i)
		BCA.write_flat_chain(output_chains[i], &cv);
	chaq::free_chaq_vecs2(cv);
	BCA.Close();
	CloseStdioFile(fOut);
	}
