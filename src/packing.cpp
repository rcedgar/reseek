#include "myutils.h"
#include "flat_chain.h"
#include "chaq.h"

void cmd_packing()
	{
	asserta(optset_maxval);
	vector<flat_chain_t *>chains;
	read_flat_chains(g_Arg1, chains);
	const uint nrchains = SIZE(chains);
	const uint M = 32;
	const int w = 16;
	vector<vector<float> > dists(w+1);
	uint maxsid = dist2sid(float(opt(maxval)));
	const uint MAXN = 99;
	vector<uint> counts(MAXN+1);
	uint maxn = 0;
	for (uint chain_idx = 0; chain_idx < nrchains; ++chain_idx)
		{
		const flat_chain_t* chain = chains[chain_idx];
		const uint L = chain->get_length();
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain->m_xyz->m_data, L, M, distmx);
		uint16_t *values = myalloc(uint16_t, L);
		chaq::get_packing_values(distmx, M, L, maxsid, true, true, values);
		for (uint i = 0; i < L; ++i)
			{
			uint n = values[i];
			maxn = max(n, maxn);
			if (n < MAXN)
				counts[n] += 1;
			}
		myfree(values);
		}

	for (uint i = 0; i <= maxn; ++i)
		ProgressLog("[%3u]  %7u\n", i, counts[i]);
	}
