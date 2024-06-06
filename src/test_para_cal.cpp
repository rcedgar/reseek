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
	Params.SetFromCmdLine(true);

	DSS D;
	D.m_Params = &Params;

	vector<vector<byte> > Profile1;
	vector<vector<byte> > Profile2;

	vector<byte> ComboLetters1;
	vector<byte> ComboLetters2;

	D.Init(Chain1);
	D.GetComboLetters(ComboLetters1);

	D.Init(Chain2);
	D.GetComboLetters(ComboLetters2);

	DSSAligner DA;
	DA.m_Params = &Params;
	DA.SetQuery(Chain1, 0, 0, &ComboLetters1);
	DA.SetTarget(Chain2, 0, 0, &ComboLetters2);

	uint Lo1, Lo2;
	string Path;
	float Score = DA.AlignComboQP_Para_Path(Lo1, Lo2, Path);
	int Score2 = DA.GetComboDPScorePathInt(ComboLetters1, ComboLetters2,
	  Lo1, Lo2, Path);
	Log("%.1f (%d) %u, %u %s\n", Score, Score2, Lo1, Lo2, Path.c_str());
	}
