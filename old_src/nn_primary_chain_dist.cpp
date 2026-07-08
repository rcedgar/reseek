#if 0
#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"

static const int s_MinDist = 12;
static const int N = 200;
static const int K = 2*N + 1;

static void IncCount(uint Pos, uint NN, 
	vector<uint> &NNDistToCount)
	{
	if (NN == UINT_MAX)
		return;
	int d = int(NN) - int(Pos);
	if (abs(d) < s_MinDist)
		return;
	int k = d + N;
	if (k >= 0 && k < NNDistToCount.size())
		NNDistToCount[k] += 1;
	}

void cmd_nn_primary_chain_dist()
	{
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	FILE *f = CreateStdioFile(opt(output));
	const uint ChainCount = SIZE(Chains);
	DSS D;
	vector<uint> NENDistToCount(K);
	vector<uint> RENDistToCount(K);
	vector<uint> PENDistToCount(K);
	vector<uint> MENDistToCount(K);
	fprintf(f, "d\tNEN\tREN\tPEN\tMEN\n");
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *Chains[ChainIdx];
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint NEN = D.GetNEN(Pos);
			uint REN = D.GetREN(Pos);
			uint PEN = D.GetPEN(Pos);
			uint MEN = D.GetMEN(Pos);
			IncCount(Pos, NEN, NENDistToCount);
			IncCount(Pos, REN, RENDistToCount);
			IncCount(Pos, PEN, PENDistToCount);
			IncCount(Pos, MEN, MENDistToCount);
			}
		}

	for (int d = -N; d <= N; ++d)
		{
		int k = d + N;
		fprintf(f, "%d", d);
		fprintf(f, "\t%u", NENDistToCount[k]);
		fprintf(f, "\t%u", RENDistToCount[k]);
		fprintf(f, "\t%u", PENDistToCount[k]);
		fprintf(f, "\t%u", MENDistToCount[k]);
		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	}
#endif // 0