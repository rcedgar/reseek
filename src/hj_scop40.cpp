#include "myutils.h"
#include "scop40bench.h"
#include "peaker.h"

static DSSParams s_Params;
static SCOP40Bench s_SB;

static void SetFeatureWeight(FEATURE F, float w)
	{
	const uint n = SIZE(s_Params.m_Features);
	for (uint i = 0; i < n; ++i)
		{
		if (s_Params.m_Features[i] == F)
			{
			s_Params.m_Weights[i] = w;
			return;
			}
		}
	asserta(false);
	}

static double EvalScop40(const Peaker &P, const vector<double> &xv)
	{
	vector<string> VarNames;
	P.GetVarNames(VarNames);
	const uint VarCount = P.GetVarCount();
	asserta(SIZE(xv) == VarCount);
	asserta(SIZE(VarNames) == VarCount);

	const float MaxFPR = 0.005f;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const string &VarName = VarNames[VarIdx];
		if (VarName == "gapopen")	{ s_Params.m_GapOpen = -float(xv[VarIdx]); continue; }
		if (VarName == "gapext")	{ s_Params.m_GapExt = -float(xv[VarIdx]); continue; }
		if (VarName == "dpw")		{ s_Params.m_dpw = -float(xv[VarIdx]); continue; }
		if (VarName == "lddtw")		{ s_Params.m_lddtw = -float(xv[VarIdx]); continue; }
		if (VarName == "ladd")		{ s_Params.m_ladd = -float(xv[VarIdx]); continue; }
		if (VarName == "revtsv")	{ s_Params.m_revtsw = -float(xv[VarIdx]); continue; }
#define F(x)	\
		if (VarName == #x)			{ SetFeatureWeight(FEATURE_##x, float(xv[VarIdx]));	continue; }
#include "featurelist.h"
		}
	s_Params.NormalizeWeights();
	//s_Params.LoadFeatures();

	s_SB.ClearHits();
	s_SB.m_Level = "sf";
	s_SB.m_ScoresAreEvalues = false;
	s_SB.RunPrealigned(opt(input2));
	s_SB.SetStats(MaxFPR);
	return s_SB.m_Area;
	}

void cmd_hj_scop40()
	{
	optset_fast = true;
	optused_fast = true;
	asserta(optset_hjspec);
	s_Params.SetDSSParams(DM_UseCommandLineOption, SCOP40_DBSIZE);
	s_SB.m_Params = &s_Params;
	s_SB.LoadDB(g_Arg1);
	s_Params.m_DBSize = (float) s_SB.GetDBSize();
	s_SB.m_QuerySelf = true;
	s_SB.m_ScoresAreEvalues = true;
	s_SB.Setup();

	vector<string> SpecLines;
	ReadLinesFromFile(opt(hjspec), SpecLines);

	Peaker P;
	asserta(optset_output);
	P.m_fTsv = CreateStdioFile(opt(output));
	setbuf(P.m_fTsv, 0);
	P.Init(SpecLines, EvalScop40);
	P.Run();
	}
