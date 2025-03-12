#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"
#include "kmeans.h"
#include "trainer.h"

static vector<int> ivalues;
static vector<int> jvalues;

// All-vs-all distances in range -r .. +r
//  excluding adjacent pairs
static void Setijs(int r)
	{
	for (int i = -r; i <= r; ++i)
		{
		for (int j = i+1; j <= r; ++j)
			{
			int minij = min(i, j);
			int maxij = max(i, j);
			if (maxij - minij == 1)
				continue;
			ivalues.push_back(minij);
			jvalues.push_back(maxij);
			}
		}
	}

static void Getv(const PDBChain &Chain, uint Pos, vector<float> &v)
	{
	v.clear();
	const uint L = Chain.GetSeqLength();
	const uint M = SIZE(ivalues);
	asserta(SIZE(jvalues) == M);
	int IntPos = Pos;
	for (uint m = 0; m < M; ++m)
		{
		int i = ivalues[m];
		int j = jvalues[m];
		float d = Chain.GetDist(IntPos+i, IntPos+j);
		v.push_back(d);
		}
	}

static void Getvs(const vector<PDBChain *> &Chains, uint N, uint r,
				  vector<vector<float> > &vs)
	{
	vs.clear();
	const uint ChainCount = SIZE(Chains);
	N = min(N, ChainCount);
	vs.reserve(N);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		const PDBChain &Chain = *Chains[ChainIndex];
		const string &Label = Chain.m_Label;

		const uint L = Chain.GetSeqLength();
		for (int Pos = r; Pos + r < int(L); ++Pos)
			{
			vector<float> v;
			Getv(Chain, Pos, v);
			if (SIZE(vs) == N)
				return;
			vs.push_back(v);
			}
		}
	}

void cmd_sscluster2()
	{
// Max number of positions to consider, for efficiency
	const uint N = optset_n ? opt(n) : 100000;
	asserta(optset_k);
	uint K = opt(k);
	const int r = (optset_r ? opt(r) : 2);
	ProgressLog("K=%u, N=%u, r=%d\n", K, N, r);

	Trainer Tr;
	Tr.Init(g_Arg1, opt(db));
	vector<PDBChain *> Chains = Tr.m_Chains;
	const uint ChainCount = SIZE(Chains);

	Setijs(r);

	vector<vector<float> > vs;
	Getvs(Chains, N, r, vs);

	KMeans KM;
	const uint MaxIters = 1000;
	KM.Run(vs, K, MaxIters);
	}
