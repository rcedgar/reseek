#include "myutils.h"
#include "dss.h"
#include "chaq.h"
#include "flat_chain.h"
#include "pdbchain.h"
#include "fast_dist_mx.h"
#include "pdbfilescanner.h"
#include "flat_chain_reader.h"

static void test_dist_mx(DSS &D, chaq &c)
	{
	const flat_chain *chain = c.m_chain;
	const PDBChain &Chain = *D.m_Chain;
	const chaindistmx_t *distmx = c.get_distmx();
	const uint16_t *distmx_ics = distmx->m_data;
	uint L = D.GetSeqLength();
	asserta(chain->get_length() == L);
	const int Li = L;
	const uint band_size = band_K(L);
	uint band_counter = 0;
	uint same = 0;
	uint diff1 = 0;
	uint diffgt1 = 0;
	for (int i = 0; i < Li; ++i)
		{
		for (int j = 0; j < i; ++j)
			{
			if (i-j >= dist_mx_band_width)
				continue;
			++band_counter;
			float d = Chain.GetDist(uint(i), uint(j));
			uint16_t dIC = Chain.CoordToIC(d);
			uint k = band_ij_to_k(i, j);
			asserta(k < band_size);
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
		same, GetPct(same, band_counter),
		diff1, GetPct(diff1, band_counter),
		diffgt1, GetPct(diffgt1, band_counter),
		chain->m_label.c_str());
//	asserta(band_size == band_counter);
	}

static void test_nn(DSS &D, chaq &c)
	{
	const flat_chain *chain = c.m_chain;
	const PDBChain &Chain = *D.m_Chain;
	c.get_distmx();

	const uint L = Chain.GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		uint PEN = D.GetPEN(i);
		uint MEN = D.GetMEN(i);
		uint pen = c.get_pen(i);
		uint men = c.get_men(i);
		Log("%4u | %4u  %4u | %4u  %4u\n", i, PEN, pen, MEN, men);
		}
	}

void cmd_test_flat()
	{
	vector<PDBChain *> Chains;
	vector<flat_chain *> chains;
	ReadChains(g_Arg1, Chains);
	const uint ChainCount = SIZE(Chains);

	DSS D;
	chaq c;

	PDBFileScanner FS;
	FS.Open(g_Arg1);
	flat_chain_reader CR;
	CR.Open(FS);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Processing");
		flat_chain *chain = CR.GetNext();
		asserta(chain);
		const PDBChain &Chain = *Chains[ChainIdx];
		asserta(chain->m_label == Chain.m_Label);

		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		if (L < 8)
			continue;

		c.init(chain);
		test_nn(D, c);
		//test_dist_mx(D, c);
		_chkmem();

		delete chain;
		}
	log_flat_stats();
	}