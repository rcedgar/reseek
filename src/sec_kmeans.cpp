#include "myutils.h"
#include "flat_chain.h"
#include "logodds.h"
#include "trainer.h"
#include "sort.h"
#include "sec_cluster.h"
#include "chaq.h"
#include "quarts.h"

void cmd_sec_kmeans()
	{
	vector<flat_chain *> chains;
	read_flat_chains(g_Arg1, chains);
	const uint nrchains = SIZE(chains);

/***
     _________________________________________
	|     |     |     | [7] |     | o-o |  o  |
	|-3,+3|-2,+3|-1,+3| 0,+3|+1,+3|+2,+3|+3,+3|    o = zero
	|_____|_____|_____|_____|_____|_____|_____|    o-o = 3.81 A (adjacent)
	|     | [2] | [4] | [5] | o-o |  o  | o-o |
	|-3,+2|-2,+2|-1,+2| 0,+2|+1,+2|+2,+2|+3,+2|
	|_____|_____|_____|_____|_____|_____|_____|
	|     | [1] | [3] | o-o |  o  | o-o |     |
	|-3,+1|-2,+1|-1,+1| 0,+1|+1,+1|+2,+1|+3,+1|
	|_____|_____|_____|_____|_____|_____|_____|
	| [8] | [0] | o-o |  o  | o-o |     |     |
	|-3,0 |-2,0 |-1,0 | 0,0 |+1,0 |+2,0 |+3,0 |
	|_____|_____|_____|_____|_____|_____|_____|
	|     | o-o |  o  | o-o |     |     |     |
	|-3,-1|-2,-1|-1,-1| 0,-1|+1,-1|+2,-1|+3,-1|
	|_____|_____|_____|_____|_____|_____|_____|
	| o-o |  o  | o-o |     |     |     |     |
	|-3,-2|-2,-2|-1,-2| 0,-2|+1,-2|+2,-2|+3,-2|
	|_____|_____|_____|_____|_____|_____|_____|
	|  o  | o-o |     |     |     |     | [6] |
	|-3,-3|-2,-3|-1,-3| 0,-3|+1,-3|+2,-3|+3,-3|
	|_____|_____|_____|_____|_____|_____|_____|

***/
	//                     0   1   2   3   4   5   6   7   8
	const vector<int> off1s = { -2, -2, -2, -1, -1,  0, -3,  0, -3 };
	const vector<int> off2s = {  0,  1,  2,  1,  2,  2,  3,  3,  0 };
	const int w = 3;

	const uint K = 16;
	const uint32_t M = 32; // dist mx band width

	sec_cluster SC;
	SC.init(K, M, off1s, off2s);

	sid_t *vs = myalloc(sid_t, SC.m_N*SC.m_D);
	SC.m_cluster_idxs = myalloc(uint, SC.m_N);
	SC.set_vs(chains);
	SC.assign_random_means();
	SC.assign_clusters();
	SC.logme();

	const uint ITERS = 1000;
	for (uint iter = 0; iter < ITERS; ++iter)
		{
		SC.iter();
		ProgressLog("iter %u, zero %u, changes %u\n",
			iter, SC.m_zero_count, SC.m_nrchanges);
		if (SC.m_nrchanges == 0)
			{
			ProgressLog("Converged\n");
			break;
			}
		}
	SC.logme();
	}

static bool check_backbone(const flat_chain *chain)
	{
	const uint L = chain->get_length();
	for (uint i = 1; i < L; ++i)
		{
		float d = chain->slow_float_dist(i-1, i);
		if (d < 3.7 || d > 3.9)
			return false;
		}
	return true;
	}

void cmd_sec_variance()
	{
	asserta(optset_output);
	FILE *fOut = CreateStdioFile(opt(output));
	vector<flat_chain *> chains;
	read_flat_chains(g_Arg1, chains);
	const uint nrchains = SIZE(chains);
	chaq c;
	const uint M = 32;
	const int w = 16;
	vector<vector<float> > dists(w+1);
	for (uint chain_idx = 0; chain_idx < nrchains; ++chain_idx)
		{
		const flat_chain *chain = chains[chain_idx];
		if (!check_backbone(chain))
			continue;
		const uint L = chain->get_length();
		for (uint pos = 0; pos < L; ++pos)
		for (uint dij = 1; dij < w; ++dij)
			{
			uint pos2 = pos + dij;
			if (pos2 >= L)
				break;
			float d = chain->slow_float_dist(pos, pos2);
			dists[dij].push_back(d);
			}
		}

	fprintf(fOut, "dij");
	fprintf(fOut, "\tN");
	fprintf(fOut, "\tMin");
	fprintf(fOut, "\tLoQ");
	fprintf(fOut, "\tMed");
	fprintf(fOut, "\tHiQ");
	fprintf(fOut, "\tMax");
	fprintf(fOut, "\tAvg");
	fprintf(fOut, "\tStdDev");
	fprintf(fOut, "\n");
	for (uint dij = 1; dij < w; ++dij)
		{
		const vector<float> &ds = dists[dij];
		QuartsFloat QF;
		GetQuartsFloat(ds, QF);
		ProgressLog("dij=%2u  ", dij);
		QF.ProgressLogMe();

		fprintf(fOut, "%u\t", dij);
		QF.ToTsv(fOut);
		}
	CloseStdioFile(fOut);
	}
