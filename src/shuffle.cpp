#if 0 // @@DELETE
#include "myutils.h"
#include "pdbchain.h"
#include "bcadata.h"

void cmd_shuffle()
	{
	asserta(optset_bca);
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	const uint N = SIZE(Chains);

	BCAData BCA;
	BCA.Create(opt(bca));
	vector<uint> v;
	for (uint i = 0; i < N; ++i)
		v.push_back(i);
	Shuffle(v);

	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Shuffling");
		PDBChain &Chain = *Chains[v[i]];
		BCA.WriteChain(Chain);
		}
	BCA.Close();
	}
#endif