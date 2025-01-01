#include "myutils.h"
#include "dssaligner.h"
#include "pdbchain.h"

void cmd_mukmerfilter()
	{
	opt_fast = true;
	optset_fast = true;
	DSSParams Params;
	Params.SetFromCmdLine(10000);
	asserta(optset_output);
	FILE *fOut = CreateStdioFile(opt_output);
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	const uint N = SIZE(Chains);
	vector<vector<byte> > MuLettersVec(N);
	vector<vector<uint> > MuKmersVec(N);
	DSS D;
	D.SetParams(Params);
	vector<string> Labels;
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Mu k-mers");
		const PDBChain &Chain = *Chains[i];
		Labels.push_back(Chain.m_Label);
		D.Init(Chain);
		D.GetMuLetters(MuLettersVec[i]);
		D.GetMuKmers(MuLettersVec[i], MuKmersVec[i]);
		}

	const uint PairCount = (N*(N-1))/2;
	uint PairIndex = 0;
	MuKmerFilter MKF;
	MKF.SetParams(Params);
	for (uint i = 0; i < N; ++i)
		{
		const string &Labeli = Labels[i];
		MKF.ResetQ();
		MKF.SetQ(&MuLettersVec[i], &MuKmersVec[i]);
		for (uint j = 0; j < i; ++j)
			{
			ProgressStep(PairIndex++, PairCount, "Filtering");
			const string &Labelj = Labels[j];
			MKF.Align(MuLettersVec[j], MuKmersVec[j]);
			if (MKF.m_BestHSPScore > 0)
				fprintf(fOut, "%s\t%s\t%d\t%d\n", Labeli.c_str(), Labelj.c_str(),
						MKF.m_BestHSPScore, MKF.m_BestChainScore);
			}
		}
	CloseStdioFile(fOut);
	}
