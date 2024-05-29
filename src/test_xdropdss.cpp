#include "myutils.h"
#include "xdpmem.h"
#include "scop40bench.h"
#include "cigar.h"

void PrettyAln(FILE *f, const PDBChain &A, const PDBChain &B,
  uint LoA, uint LoB, const string &Path, float Evalue);

static const vector<vector<byte> > *ptrProfileA;
static const vector<vector<byte> > *ptrProfileB;
static DSSAligner DA;

static float SubFn(void *UserData, uint PosA, uint PosB)
	{
	float Score = DA.GetScorePosPair(*ptrProfileA, *ptrProfileB, PosA, PosB);
	return Score;
	}

void cmd_test_xdropdss()
	{
	const string &CalFN = g_Arg1;
	SCOP40Bench SB;
	DSSParams Params;
	Params.SetDefaults();
	SB.Setup(Params);
	SB.ReadChains(CalFN);

	DSS DKF;
	DKF.m_Params = &Params;
	DA.m_Params = &Params;

	const uint i = 0;
	const uint j = 1;
	asserta(i < SIZE(SB.m_Chains));
	asserta(j < SIZE(SB.m_Chains));
	const PDBChain &ChainA = *SB.m_Chains[i];
	const PDBChain &ChainB = *SB.m_Chains[j];
	const string &LabelA = ChainA.m_Label;
	const string &LabelB = ChainB.m_Label;
	const uint LA = ChainA.GetSeqLength();
	const uint LB = ChainB.GetSeqLength();
	const vector<vector<byte> > &ProfileA = SB.m_Profiles[i];
	const vector<vector<byte> > &ProfileB = SB.m_Profiles[j];
	ptrProfileA = &ProfileA;
	ptrProfileB = &ProfileB;

	DA.SetQuery(ChainA, ProfileA, 0, 0);
	DA.SetTarget(ChainB, ProfileB, 0, 0);
	DA.Align_NoAccel();
	string CIGAR;
	PathToCIGAR(DA.m_PathAB.c_str(), CIGAR);
	float SWPathScore = DA.GetDPScorePath(ProfileA, ProfileB,
	  DA.m_LoA, DA.m_LoB, DA.m_PathAB);

	ProgressLog("\n");
	ProgressLog(">%s, %s\n", LabelA.c_str(), LabelB.c_str());
	ProgressLog("SWpathscore %.3g\n", SWPathScore);
	ProgressLog("SW E-value  %.3g\n", DA.m_EvalueAB);
	ProgressLog("SW CIGAR    (%u,%u) %s\n",
	  DA.m_LoA, DA.m_LoB, CIGAR.c_str());

	vector<uint> KmersA;
	vector<uint> KmersB;
	DKF.Init(ChainA);
	DKF.GetComboKmers(KmersA);
	DKF.Init(ChainB);
	DKF.GetComboKmers(KmersB);

	vector<uint> DiagPosAs;
	vector<uint> DiagPosBs;
	vector<uint> DiagLengths;
	DA.GetDiagSeedPairs(ProfileA, ProfileB, KmersA, KmersB,
	  DiagPosAs, DiagPosBs, DiagLengths);

	vector<int> Diags;
	const uint n = SIZE(DiagLengths);
	//if (n == 0)
	//	{
	//	ProgressLog("No diags\n");
	//	return;
	//	}

	asserta(SIZE(DiagPosAs) == n);
	asserta(SIZE(DiagPosBs) == n);
	for (uint i = 0; i < n; ++i)
		{
		uint Diag = int(DiagPosAs[i]) - int(DiagPosBs[i]);
		Diags.push_back(Diag);
		}

	float X = 2;
	float Open = -1.38f;
	float Ext = -0.065f;
	XDPMem Mem;
	string XPath;
	string XCIGAR;
	string XCIGARB;
	uint LoA = DiagPosAs[15];
	uint LoB = DiagPosBs[15];
	float XScore = XDropFwd(Mem, X, Open, Ext, SubFn, 0, LoA, LA, LoB, LB, XPath);
	float XPathScore = DA.GetDPScorePath(ProfileA, ProfileB, LoA, LoB, XPath);

	string XPathBwd;
	float XScoreBwd = XDropBwd(Mem, X, Open, Ext, SubFn, 0, LoA, LA, LoB, LB, XPathBwd);
	float XPathScoreBwd = DA.GetDPScorePath(ProfileA, ProfileB, LoA, LoB, XPathBwd);
	uint LoLoA = LoA + 1;
	uint LoLoB = LoB + 1;
	uint ColCountB = SIZE(XPathBwd);
	for (uint Col = 0; Col < ColCountB; ++Col)
		{
		char c = XPathBwd[Col];
		if (c == 'M' || c == 'D')
			{
			asserta(LoLoA > 0);
			--LoLoA;
			}
		if (c == 'M' || c == 'I')
			{
			asserta(LoLoB > 0);
			--LoLoB;
			}
		}

	PathToCIGAR(XPath.c_str(), XCIGAR);
	PathToCIGAR(XPathBwd.c_str(), XCIGARB);
	ProgressLog("X score     %.3g\n", XScore);
	ProgressLog("X scorebwd  %.3g\n", XScoreBwd);
	ProgressLog("Xpathscore  %.3g\n", XPathScore);
	ProgressLog("XpathscoreB %.3g\n", XPathScoreBwd);
	ProgressLog("XCIGAR      (%u,%u) %s\n", LoA, LoB, XCIGAR.c_str());
	ProgressLog("XCIGARB     (%u,%u) %s\n", LoLoA, LoLoB, XCIGARB.c_str());

	PrettyAln(g_fLog, ChainA, ChainB, DA.m_LoA, DA.m_LoB, DA.m_PathAB, -1);
	PrettyAln(g_fLog, ChainA, ChainB, LoA, LoB, XPath, -1);
	}
