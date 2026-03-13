#include "myutils.h"
#include "flat_chain.h"

static void get_codeseq(
	const flat_chain_t *chain,
	const string &feature,
	vector<uint8_t> &codeseq)
	{
	codeseq.clear();
	const uint L = chain->get_length();
	codeseq.resize(L);

	Die("TODO");
	}

void cmd_flat_quantize()
	{
	const string &chainfn = g_Arg1;
	vector<vector<uint8_t> > codeseqs;
	asserta(optset_alpha_size);
	const uint alpha_size = opt(alpha_size);

	asserta(optset_feature);
	const string feature = opt(feature);

	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	uint nchain = SIZE(chains);

	vector<uint8_t> codes;
	uint nbad = 0;
	for (uint i = 0; i < nchain; ++i)
		{
		vector<uint8_t> codeseq;
		const flat_chain_t *chain = chains[i];
		get_codeseq(chain, feature, codeseq);
		for (uint j = 0; j < SIZE(codeseq); ++j)
			{
			uint8_t code = codeseq[j];
			if (code < alpha_size)
				codes.push_back(code);
			else
				++nbad;
			}
		}
	ProgressLog("%u good, %u bad codes\n", SIZE(codes), nbad);
	sort(codes.begin(), codes.end());

	const uint N = SIZE(codes);
	vector<uint> ts;
	for (uint i = 0; i + 1 < alpha_size; ++i)
		{
		uint k = (i*N)/alpha_size;
		ts.push_back(codes[k]);
		}

	for (uint i = 0; i + 1 < alpha_size; ++i)
		ProgressLog("ts[%2u] = %u\n", i, ts[i]);
	}
