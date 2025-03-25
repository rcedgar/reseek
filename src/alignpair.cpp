#include "myutils.h"
#include "chainreader2.h"
#include "dssaligner.h"
#include "pdbchain.h"
#include "xyz.h"
#include "dss.h"

float GetSelfRevScore(DSSAligner &DA, const DSSParams &Params,
					  const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers)
	{
	if (opt(selfrev0))
		return 0;

	DSS D;
	PDBChain RevChain;
	Chain.GetReverse(RevChain);
	vector<vector<byte> > RevProfile;
	D.Init(RevChain, Params);
	D.GetProfile(RevProfile);

	DA.SetQuery(Chain, &Profile, ptrMuLetters, ptrMuKmers, FLT_MAX);

	DA.SetTarget(RevChain, &RevProfile, ptrMuLetters, ptrMuKmers, FLT_MAX);
	DA.AlignQueryTarget();
	return DA.m_AlnFwdScore;
	}

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

static void XformLine(const float t[3],
  const float u[3][3], string &Line)
	{
	float x, y, z;
	PDBChain::GetXYZFromATOMLine(Line, x, y, z);

	float Pt[3];
	float XPt[3];

	Pt[0] = x;
	Pt[1] = y;
	Pt[2] = z;
	transform(t, u, Pt, XPt);

	x = (float) XPt[0];
	y = (float) XPt[1];
	z = (float) XPt[2];
	PDBChain::SetXYZInATOMLine(Line, x, y, z, Line);
	}

static void XformLines(const float t[3],
  const float u[3][3], vector<string> &Lines)
	{
	const uint N = SIZE(Lines);
	for (uint i = 0; i < N; ++i)
		{
		string &Line = Lines[i];
		if (PDBChain::IsATOMLine(Line))
			XformLine(t, u, Line);
		}
	}

static float AlignPair1(const DSSParams &Params, DSS &D, DSSAligner &DA,
  const PDBChain *ChainQ, const PDBChain *ChainT, bool DoOutput)
	{
	vector<vector<byte> > ProfileQ;
	vector<vector<byte> > ProfileT;

	vector<byte> MuLettersQ;
	vector<uint> MuKmersQ;

	float BestScore = 0;
	uint BestChainIndexQ = UINT_MAX;
	uint BestChainIndexT = UINT_MAX;
	D.Init(*ChainQ, Params);
	D.GetProfile(ProfileQ);
	if (Params.m_Omega > 0)
		D.GetMuLetters(MuLettersQ);

	float SelfRevScoreQ = GetSelfRevScore(DA, Params, *ChainQ, ProfileQ, 0, 0);

	vector<byte> MuLettersT;
	vector<uint> MuKmersT;

	D.Init(*ChainT, Params);
	D.GetProfile(ProfileT);

	if (Params.m_Omega > 0)
		D.GetMuLetters(MuLettersT);

	float SelfRevScoreT = GetSelfRevScore(DA, Params, *ChainT, ProfileT, 0, 0);
	
	DA.SetQuery(*ChainQ, &ProfileQ, &MuLettersQ, &MuKmersQ, SelfRevScoreQ);
	DA.SetTarget(*ChainT, &ProfileT, &MuLettersT, &MuKmersT, SelfRevScoreT);
	float Score = FLT_MAX;
	if (opt(global))
		{
		DA.AlignQueryTarget_Global();
		Score = DA.m_GlobalScore;
		}
	else
		{
		DA.AlignQueryTarget();
		Score = DA.m_AlnFwdScore;
		}
	if (DoOutput)
		{
		if (optset_aln)
			{
			FILE *f = CreateStdioFile(opt(aln));
			DA.ToAln(f, true);
			CloseStdioFile(f);
			}

		float t[3];
		float u[3][3];
		DA.GetKabsch(t, u, true);

		vector<string> LinesQ = ChainQ->m_Lines;
		XformLines(t, u, LinesQ);

		if (optset_output)
			{
			FILE *f = CreateStdioFile(opt(output));
			for (uint i = 0; i < SIZE(LinesQ); ++i)
				fprintf(f, "%s\n", LinesQ[i].c_str());
			CloseStdioFile(f);
			}
		}
	return Score;
	}

void cmd_alignpair()
	{
	if (!optset_input2)
		Die("Must specify -input2");

	const string &QFN = g_Arg1;
	const string &TFN = opt(input2);

	vector<PDBChain *> ChainsQ;
	vector<PDBChain *> ChainsT;
	ReadChains_SaveLines(QFN, ChainsQ);
	ReadChains_SaveLines(TFN, ChainsT);

	optset_sensitive = true;
	opt(sensitive) = true;
	DSSParams Params;
	Params.SetDSSParams(DM_AlwaysSensitive, SCOP40_DBSIZE);
	Params.m_UsePara = false;
	Params.m_Omega = 0;

	DSSAligner DA;
	DA.SetParams(Params);

	DSS D;
	const uint ChainCountQ = SIZE(ChainsQ);
	const uint ChainCountT = SIZE(ChainsT);
	if (ChainCountQ == 0) Die("No chains found in %s", QFN.c_str());
	if (ChainCountT == 0) Die("No chains found in %s", TFN.c_str());

	vector<vector<byte> > ProfileQ;
	vector<vector<byte> > ProfileT;

	vector<byte> MuLettersQ;
	vector<uint> KmersQ;

	float BestScore = -9999;
	uint BestChainIndexQ = UINT_MAX;
	uint BestChainIndexT = UINT_MAX;
	for (uint ChainIndexQ = 0; ChainIndexQ < ChainCountQ; ++ChainIndexQ)
		{
		PDBChain *ChainQ = ChainsQ[ChainIndexQ];
		for (uint ChainIndexT = 0; ChainIndexT < ChainCountT; ++ChainIndexT)
			{
			PDBChain *ChainT = ChainsT[ChainIndexT];
			float Score = AlignPair1(Params, D, DA, ChainQ, ChainT, false);
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
	PDBChain *ChainT = ChainsT[BestChainIndexT];
	AlignPair1(Params, D, DA, ChainQ, ChainT, true);
	}
