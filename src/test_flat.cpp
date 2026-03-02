#include "myutils.h"
#include "dss.h"
#include "chaq.h"
#include "flat_chain.h"
#include "pdbchain.h"
#include "fast_dist_mx2.h"

void cmd_test_flat()
	{
	vector<PDBChain *> Chains;
	vector<flat_chain *> chains;
	ReadChains(g_Arg1, Chains);
	read_flat_chains(g_Arg1, chains);
	const uint ChainCount = SIZE(Chains);
	asserta(SIZE(chains) == ChainCount);

	DSS D;
	chaq c;

	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Processing");
		const PDBChain &Chain = *Chains[ChainIdx];
		const flat_chain *chain = chains[ChainIdx];
		D.Init(Chain);
		c.init(chain);
		const chaindistmx_t *distmx = c.get_distmx();
		const uint16_t *distmx_ics = distmx->m_data;

		const uint L = Chain.GetSeqLength();
		asserta(chain->get_length() == L);
		const int Li = L;
		const uint band_size = band_K(L, M);
		uint band_counter = 0;
		uint same = 0;
		uint diff1 = 0;
		uint diffgt1 = 0;
		for (int i = 0; i < Li; ++i)
			{
			for (int j = 0; j < i; ++j)
				{
				if (i-j >= M)
					continue;
				++band_counter;
				float d = Chain.GetDist(uint(i), uint(j));
				uint16_t dIC = Chain.CoordToIC(d);
				uint k = band_ij_to_k(i, j, L, M);
				uint16_t dIC2 = distmx_ics[k];
				int diff = int(dIC2) - int(dIC);
				if (diff == 0)
					++same;
				else if (abs(diff) == 1)
					++diff1;
				else
					++diffgt1;
				//Log("%d", i);
				//Log("\t%d", j);
				//Log("\t%u", dIC);
				//Log("\t%u", dIC2);
				//if (diff != 0)
				//	Log("\t%d", diff);
				//Log("\n");
				}
			}
		Log("same=%7u (%5.1f%%), diff1=%7u (%5.1f%%), diffgt1=%7u (%5.1f%%) %s\n",
			same, GetPct(same, band_size),
			diff1, GetPct(diff1, band_size),
			diffgt1, GetPct(diffgt1, band_size),
			chain->m_label.c_str());
		asserta(band_size == band_counter);
		}
	_chkmem();
	log_flat_stats();
	}