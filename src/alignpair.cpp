#include "myutils.h"
#include "chainreader2.h"
#include "dssaligner.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "alncounts.h"

float GetSelfRevScore(DSSAligner &DA, DSS &D, const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers);

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

static float AlignPair1(DSS &D, DSSAligner &DA,
  const PDBChain *ChainQ, const PDBChain *ChainT, bool DoOutput)
	{
	vector<vector<byte> > ProfileQ;
	vector<vector<byte> > ProfileT;

	vector<byte> MuLettersQ;
	vector<uint> MuKmersQ;

	float BestScore = 0;
	uint BestChainIndexQ = UINT_MAX;
	uint BestChainIndexT = UINT_MAX;
	D.Init(*ChainQ);
	D.GetProfile(ProfileQ);
	int Omega = DSSParams::GetOmega();
	if (Omega > 0)
		D.GetMuLetters(MuLettersQ);

	float SelfRevScoreQ = GetSelfRevScore(DA, D, *ChainQ, ProfileQ, 0, 0);

	vector<byte> MuLettersT;
	vector<uint> MuKmersT;

	D.Init(*ChainT);
	D.GetProfile(ProfileT);

	if (Omega > 0)
		D.GetMuLetters(MuLettersT);

	float SelfRevScoreT = GetSelfRevScore(DA, D, *ChainT, ProfileT, 0, 0);
	
	DA.SetQuery(*ChainQ, &ProfileQ, &MuLettersQ, &MuKmersQ, SelfRevScoreQ);
	DA.SetTarget(*ChainT, &ProfileT, &MuLettersT, &MuKmersT, SelfRevScoreT);
	DA.AlignQueryTarget();
	float Score = DA.m_AlnFwdScore;
	if (DoOutput)
		{
		if (optset_aln)
			{
			FILE *f = CreateStdioFile(opt(aln));
			DA.ToAln(f, true);
			CloseStdioFile(f);
			}

		double t[3];
		double u[3][3];
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
		if (optset_output2)
			{
			FILE *f = CreateStdioFile(opt(output2));
			for (uint i = 0; i < SIZE(LinesQ); ++i)
				{
				string Line = LinesQ[i];
				Line[21] = '1';
				fprintf(f, "%s\n", Line.c_str());
				}

			const vector<string> &LinesT = ChainT->m_Lines;
			for (uint i = 0; i < SIZE(LinesT); ++i)
				{
				string Line = LinesT[i];
				Line[21] = '2';
				fprintf(f, "%s\n", Line.c_str());
				}
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
	DSSParams::Init(DM_AlwaysSensitive);
	DSSParams::m_Omega8 = 0;
	DSSParams::m_Omega16 = 0;

	DSSAligner DA;
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
			float Score = AlignPair1(D, DA, ChainQ, ChainT, false);
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
	AlignPair1(D, DA, ChainQ, ChainT, true);
	}
