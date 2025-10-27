#include "myutils.h"
#include "statsig.h"
#include "scop40bench.h"
#include "peaker.h"

static SCOP40Bench *s_SB;
static Peaker *s_Peaker;

static double EvalArea(const vector<double> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);
	string VarsStr;
	s_Peaker->VarsToStr(xv, VarsStr, ";");
	DSSParams::SetTunableParamsFromStr(VarsStr);
	s_SB->ClearHitsAndResults();
	s_SB->RunSelf(false);
	s_SB->m_Level = "sf";
	s_SB->SetStats(0.005f);
	s_SB->WriteSummary();
	return s_SB->m_Area;
	}

void cmd_evalarea()
	{
	const string &DBFN = g_Arg1;

	DSSParams::Init(DM_UseCommandLineOption);
	string ParamStr;
	DSSParams::GetParamStr(ParamStr);
	ProgressLog("ParamStr:%s\n", ParamStr.c_str());

	s_SB = new SCOP40Bench;
	s_SB->LoadDB(DBFN);
	StatSig::Init(s_SB->GetDBSize());
	s_SB->Setup();
	s_SB->m_QuerySelf = true;
	s_SB->RunSelf(false);
	s_SB->m_Level = "sf";
	s_SB->SetStats(0.005f);
	s_SB->WriteSummary();
	ProgressLog("Area=%.4g\n", s_SB->m_Area);
	}

void cmd_hjmega()
	{
	const string SpecFN = g_Arg1;
	asserta(optset_db);
	const string &DBFN = opt(db);

	DSSParams::Init(DM_UseCommandLineOption);
	string ParamStr;
	DSSParams::GetParamStr(ParamStr);
	ProgressLog("ParamStr:%s\n", ParamStr.c_str());

	s_SB = new SCOP40Bench;
	s_SB->LoadDB(DBFN);
	StatSig::Init(s_SB->GetDBSize());
	s_SB->Setup();
	s_SB->m_QuerySelf = true;

	Log("SpecFN=%s\n", SpecFN.c_str());
	vector<string> SpecLines;
	ReadLinesFromFile(SpecFN, SpecLines);
	for (uint i = 0; i < SIZE(SpecLines); ++i)
		Log("%s\n", SpecLines[i].c_str());

	Peaker P;
	if (optset_output2)
		P.m_fTsv = CreateStdioFile(opt(output2));
	s_Peaker = &P;
	P.Init(SpecLines, EvalArea);
	//P.LogLatinBins();
	if (opt(latinloop))
		P.LatinLoop();
	else if (optset_nested_latin)
		P.RunNestedLatin(opt(nested_latin));
	P.Run();
	CloseStdioFile(P.m_fTsv);
	}
