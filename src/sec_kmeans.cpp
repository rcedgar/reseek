#include "myutils.h"
#include "flat_chain.h"
#include "logodds.h"
#include "trainer.h"
#include "sort.h"
#include "sec_cluster.h"
#include "chaq.h"

// -2,0 
// -2,1
// -2,2
// -1,1
// -1,2
// 0,2
// -3,3
// 0,3
// -3,0

void cmd_sec_kmeans()
	{
	vector<flat_chain *> chains;
	read_flat_chains(g_Arg1, chains);
	const uint nrchains = SIZE(chains);

	uint N = 0;
	for (uint chainidx = 0; chainidx < nrchains; ++chainidx)
		N += chains[chainidx]->get_length();

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
	const int offs1[] = { -2, -2, -2, -1, -1,  0, -3,  0, -3 };
	const int offs2[] = {  0,  1,  2,  1,  2,  2,  3,  3,  0 };
	const int w = 3;

	const uint32_t M = 32; // dist mx band width
	sec_cluster SC;
	SC.m_K = 16;
	SC.m_off1s = offs1;
	SC.m_off2s = offs2;
	SC.m_w = w;
	SC.m_D = 9;
	SC.m_M = M;
	SC.m_N = N;
	SC.m_means = myalloc(sid_t, SC.m_K*SC.m_D);
	sid_t *vs = myalloc(sid_t, SC.m_N*SC.m_D);
	SC.m_cluster_idxs = myalloc(uint, SC.m_N);
	SC.set_vs(chains);
	SC.assign_random_means();
	SC.assign_clusters();
	SC.logme();

	const uint ITERS = 100;
	for (uint iter = 0; iter < ITERS; ++iter)
		{
		uint zero_count = SC.calc_means();
		uint nrchanges = SC.assign_clusters();
		ProgressLog("iter %u, zero %u, changes %u\n", iter, zero_count, nrchanges);
		if (nrchanges == 0)
			{
			ProgressLog("Converged\n");
			uint nrchanges2 = SC.assign_clusters();
			asserta(nrchanges2 == 0);
			break;
			}
		}
	SC.logme();
	}
