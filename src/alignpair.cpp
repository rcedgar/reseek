"myutils.h"
#include "chainreader2.h"
#include "dssaligner.h"
#include "pdbchain.h"
#include "abcxyz.h"

static void ReadChains_SaveLines(const string &FileName,
  vector<PDBChain *> &Chains)
	{
	PDBFileScanner FS;
	FS.Open(FileName);

	ChainReader2 CR;
	CR.Open(FS);
	CR.m_SaveLines = true;
	for (;;)
		{
		PDBChain *Chain = CR.GetNext();
		if (Chain == 0)
			break;
		Chains.push_back(Chain);
		}
	}

static void XformLine(const double t[3],
  const double u[3][3], string &Line)
	{
	float x, y, z;
	PDBChain::GetXYZFromATOMLine(Line, x, y, z);

	double Pt[3];
	double XPt[3];

	Pt[0] = x;
	Pt[1] = y;
	Pt[2] = z;
	transform(t, u, Pt, XPt);

	x = (float) XPt[0];
	y = (float) XPt[1];
	z = (float) XPt[2];
	PDBChain::SetXYZInATOMLine(Line, x, y, z, Line);
	}

static void XformLines(const double t[3],
  const double u[3][3], vector<string> &Lines)
	{
	const uint N = SIZE(Lines);
	for (uint i = 0; i < N; ++i)
		{
		string &Line = Lines[i];
		if (PDBChain::IsATOMLine(Line))
			XformLine(t, u, Line);
		}
	}

void cmd_alignpair()
	{
	if (!optset_input2)
		Die("Must specify -input2");

	const string &QFN = g_Arg1;
	const string &TFN = opt_input2;

	vector<PDBChain *> ChainsQ;
	vector<PDBChain *> ChainsT;
	ReadChains_SaveLines(QFN, ChainsQ);
	ReadChains_SaveLines(TFN, ChainsT);

	optset_sensitive = true;
	opt_sensitive = true;
	DSSParams Params;
	Params.SetFromCmdLine(10000);
	Params.m_UsePara = false;
	Params.m_Omega = 0;

	DSSAligner DA;
	DA.m_Params = &Params;

	DSS D;
	D.SetParams(Params);

	const uint ChainCountQ = SIZE(ChainsQ);
	const uint ChainCountT = SIZE(ChainsT);
	if (ChainCountQ == 0) Die("No chains found in %s", QFN.c_str());
	if (ChainCountT == 0) Die("No chains found in %s", TFN.c_str());

	vector<vector<byte> > ProfileQ;
	vector<vector<byte> > ProfileT;

	float BestScore = 0;
	uint BestChainIndexQ = UINT_MAX;
	uint BestChainIndexT = UINT_MAX;
	for (uint ChainIndexQ = 0; ChainIndexQ < ChainCountQ; ++ChainIndexQ)
		{
		PDBChain *ChainQ = ChainsQ[ChainIndexQ];
		D.Init(*ChainQ);
		D.GetProfile(ProfileQ);
		DA.SetQuery(*ChainQ, &ProfileQ, 0, 0, FLT_MAX);

		for (uint ChainIndexT = 0; ChainIndexT < ChainCountT; ++ChainIndexT)
			{
			PDBChain *ChainT = ChainsT[ChainIndexT];

			D.Init(*ChainT);
			D.GetProfile(ProfileT);

			DA.SetTarget(*ChainT, &ProfileT, 0, 0, FLT_MAX);
			DA.AlignQueryTarget();
			float Score = DA.m_AlnFwdScore;
			if (Score > BestScore)
				{
				BestScore = Score;
				BestChainIndexQ = ChainIndexQ;
				BestChainIndexT = ChainIndexT;
				}
			}
		}
	if (BestScore == 0)
		Die("No alignment found");

	PDBChain *ChainQ = ChainsQ[BestChainIndexQ];
	D.Init(*ChainQ);
	D.GetProfile(ProfileQ);
	DA.SetQuery(*ChainQ, &ProfileQ, 0, 0, FLT_MAX);

	PDBChain *ChainT = ChainsT[BestChainIndexT];

	D.Init(*ChainT);
	D.GetProfile(ProfileT);

	DA.SetTarget(*ChainT, &ProfileT, 0, 0, FLT_MAX);
	DA.AlignQueryTarget();

	if (optset_aln)
		{
		FILE *f = CreateStdioFile(opt_aln);
		DA.ToAln(f, true);
		CloseStdioFile(f);
		}
	Log("%s\n", DA.m_Path.c_str());

	double t[3];
	double u[3][3];
	DA.GetKabsch(t, u, true);

	vector<string> &LinesQ = ChainQ->m_Lines;
	vector<string> &LinesT = ChainT->m_Lines;

	XformLines(t, u, LinesQ);
	//XformLines(t, u, LinesT);

	if (optset_output)
		{
		FILE *f = CreateStdioFile(opt_output);
		for (uint i = 0; i < SIZE(LinesQ); ++i)
			fprintf(f, "%s\n", LinesQ[i].c_str());
		CloseStdioFile(f);
		}

	//if (optset_output2)
	//	{
	//	FILE *f = CreateStdioFile(opt_output2);
	//	for (uint i = 0; i < SIZE(LinesT); ++i)
	//		fprintf(f, "%s\n", LinesT[i].c_str());
	//	CloseStdioFile(f);
	//	}
	}
