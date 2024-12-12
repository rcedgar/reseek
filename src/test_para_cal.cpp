#include "myutils.h"
#include "pdbchain.h"
#include "dssaligner.h"

void cmd_test_para_cal()
	{
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	asserta(SIZE(Chains) >= 2);

	const PDBChain &Chain1 = *Chains[0];
	const PDBChain &Chain2 = *Chains[1];

	DSSParams Params;
	Params.SetFromCmdLine(10000);

	DSS D;
	D.SetParams(Params);

	vector<vector<byte> > Profile1;
	vector<vector<byte> > Profile2;

	vector<byte> MuLetters1;
	vector<byte> MuLetters2;

	D.Init(Chain1);
	D.GetMuLetters(MuLetters1);

	D.Init(Chain2);
	D.GetMuLetters(MuLetters2);

	DSSAligner DA;
	DA.SetParams(Params);
	DA.SetQuery(Chain1, 0, &MuLetters1, 0, FLT_MAX);
	DA.SetTarget(Chain2, 0, &MuLetters2, 0, FLT_MAX);

	uint Lo1, Lo2;
	string Path;
	float Score = DA.AlignMuQP_Para_Path(Lo1, Lo2, Path);
	int Score2 = DA.GetMuDPScorePathInt(MuLetters1, MuLetters2,
	  Lo1, Lo2, Path);
	Log("%.1f (%d) %u, %u %s\n", Score, Score2, Lo1, Lo2, Path.c_str());
	}
