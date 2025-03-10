#include "myutils.h"
#include "scop40bench.h"
#include "peaker.h"

static DSSParams s_Params;
static SCOP40Bench s_SB;

static double EvalScop40(const vector<double> &xv)
	{
	const float MaxFPR = 0.005f;
	s_Params.m_GapOpen = -float(xv[0]);
	s_Params.m_GapExt = -float(xv[1]);

	s_SB.ClearHits();
	s_SB.m_Level = "sf";
	s_SB.RunSelf();
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
	P.m_fTsv = CreateStdioFile(opt(output));
	setbuf(P.m_fTsv, 0);
	P.Init(SpecLines, EvalScop40);
	P.Run();
	}
