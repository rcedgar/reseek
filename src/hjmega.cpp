#include "myutils.h"
#include "statsig.h"
#include "scop40bench.h"
#include "peaker.h"

static SCOP40Bench *s_SB;
static Peaker *s_Peaker;

static double EvalArea(const vector<string> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);
	string VarsStr;
	s_Peaker->xv2str(xv, VarsStr);
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

	if (optset_subsample)
		{
		uint Pct = opt(subsample);
		SCOP40Bench *Subset = new SCOP40Bench;
		s_SB->MakeSubset(*Subset, Pct);
		Subset->Setup();
		Subset->m_QuerySelf = true;
		s_SB = Subset;
		ProgressLog("Subset %u%%, %u chains\n",
			Pct, Subset->GetDBChainCount());
		}

	Log("SpecFN=%s\n", SpecFN.c_str());
	vector<string> SpecLines;
	ReadLinesFromFile(SpecFN, SpecLines);
	for (uint i = 0; i < SIZE(SpecLines); ++i)
		Log("%s\n", SpecLines[i].c_str());

	Peaker P(0, "");
	if (optset_output2)
		P.m_fTsv = CreateStdioFile(opt(output2));
	s_Peaker = &P;
	P.Init(SpecLines, EvalArea);

	bool SkipInit = P.GetGlobalBool("skipinit", false);
	uint Latin = P.GetGlobalInt("latin", 0);
	if (!SkipInit)
		P.RunInitialValues();
	if (Latin > 0)
		P.RunLatin();
	P.HJ_RunHookeJeeves();
	CloseStdioFile(P.m_fTsv);

	string BestVarStr;
	P.xv2str(P.m_Best_xv, BestVarStr);
	ProgressLog("FINAL [%.3g] %s\n", P.m_Best_y, BestVarStr.c_str());

	if (optset_input2)
		{
		s_SB = new SCOP40Bench;
		s_SB->LoadDB(opt(input2));
		StatSig::Init(s_SB->GetDBSize());
		s_SB->Setup();
		s_SB->m_QuerySelf = true;
		double Area = EvalArea(P.m_Best_xv);
		ProgressLog("FULLDB [%.3g] %s\n", Area, BestVarStr.c_str());
		}
	}
