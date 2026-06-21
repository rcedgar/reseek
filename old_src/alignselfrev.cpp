#if 0
#include "myutils.h"
#include "dssaligner.h"
#include "chainreader2.h"

void LogProfile(
	const string &label,
	vector<vector<byte> > &profile)
	{
	uint nfeat = uint(profile.size());
	uint L = uint(profile[0].size());
	Log("LogProfile(%s) L=%u nfeat=%u\n", label.c_str(), L, nfeat);
	Log("  pos  ");
	for (uint fi = 0; fi < nfeat; ++fi)
		Log(" %2u", fi);
	Log("\n");
	for (uint pos = 0; pos < L; ++pos)
		{
		Log("[%4u] ", pos);
		for (uint fi = 0; fi < nfeat; ++fi)
			Log(" %2x", profile[fi][pos]);
		Log("\n");
		}
	}

void cmd_alignselfrev()
	{
	const string &QFN = g_Arg1;
	FILE *fOut = CreateStdioFile(opt(output));

	DSSParams::Init(DM_AlwaysSensitive);
	if (optset_varstr)
		DSSParams::SetParamsFromStr(opt(varstr));

	DSSParams::m_Omega8 = 0;
	DSSParams::m_Omega16 = 0;

	DSSAligner DA;
	DSS D;

	vector<vector<byte> > Profile;
	vector<vector<byte> > RevProfile;
	ChainReader2 CR;
	CR.Open(QFN);
	uint N = 0;
	for (;;)
		{
		PDBChain *Chain = CR.GetNext();
		if (Chain == 0)
			break;
		++N;
		if (N%100 == 0)
			Progress("%u\r", N);
		PDBChain Rev = *Chain;
		Rev.Reverse();

		D.Init(*Chain);
		D.GetProfile(Profile);
		LogProfile(Chain->m_Label, Profile);

		D.Init(Rev);
		D.GetProfile(RevProfile);
		LogProfile(Chain->m_Label + ".rev", RevProfile);

		DA.SetQuery(*Chain, &Profile, 0, 0, FLT_MAX);
		DA.SetTarget(Rev, &RevProfile, 0, 0, FLT_MAX);
		DA.AlignQueryTarget();
		DA.ToTsv(fOut, true);
		}
	CloseStdioFile(fOut);
	}
#endif
