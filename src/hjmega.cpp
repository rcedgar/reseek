#include "myutils.h"
#include "statsig.h"
#include "scop40bench.h"
#include "peaker.h"

static SCOP40Bench *s_SB;
static Peaker *s_Peaker;

static double EvalArea0(const vector<string> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);
	string VarsStr;
	s_Peaker->xv2xss(xv, VarsStr);
	DSSParams::SetTunableParamsFromStr(VarsStr);
	s_SB->ClearHitsAndResults();
	s_SB->RunSelf(false);
	s_SB->m_Level = "sf";
	s_SB->SetStats(0.005f);
	s_SB->WriteSummary();
	return s_SB->m_Area0;
	}

static double EvalArea3(const vector<string> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);
	string VarsStr;
	s_Peaker->xv2xss(xv, VarsStr);
	DSSParams::SetTunableParamsFromStr(VarsStr);
	s_SB->ClearHitsAndResults();
	s_SB->RunSelf(false);
	s_SB->m_Level = "sf";
	s_SB->SetStats(0.005f);
	s_SB->WriteSummary();
	return s_SB->m_Area3;
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
	ProgressLog("Area0=%.4g, Area3=%.4g\n", s_SB->m_Area0, s_SB->m_Area3);
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
	P.Init(SpecLines, EvalArea3);

	if (!P.m_InitParams.empty())
		{
		vector<string> xv;
		P.xss2xv(P.m_InitParams, xv);
		P.Evaluate(xv, "init");
		}
	uint Latin = P.GetGlobalInt("latin", 0);
	if (Latin > 0)
		P.RunLatin();

	vector<uint> TopEvalIdxs;
	P.GetTopEvalIdxs(3, TopEvalIdxs);
	const uint n = SIZE(TopEvalIdxs);
	if (n == 0)
		Die("No evals");
	for (uint k = 0; k < n; ++k)
		{
		ProgressLog("=========================================\n");
		ProgressLog("               HJ %u/%u\n", k+1, n);
		ProgressLog("=========================================\n");

		uint EvalIdx = TopEvalIdxs[k];
		string ChildName;
		Ps(ChildName, "HJ%u/%u", k+1, n);
		Peaker *Child = P.MakeChild(ChildName);
		double y = P.m_ys[EvalIdx];
		const vector<string> &xv = P.m_xvs[EvalIdx];
		Child->AppendResult(xv, y, "HJstart");
		Child->HJ_RunHookeJeeves();
		P.AppendChildResults(*Child);
		delete Child;
		}

	ProgressLog("Eval cache hits %u\n", P.m_EvaluateCacheHits);

	string BestVarStr;
	P.xv2xss(P.m_Best_xv, BestVarStr);

	ProgressLog("\n");
	ProgressLog("CONVERGED [%.3g] %s\n", P.m_Best_y, BestVarStr.c_str());
	ProgressLog("\n");

	P.WriteFinalResults(g_fLog);

	if (optset_input2)
		{
		s_SB = new SCOP40Bench;
		s_SB->LoadDB(opt(input2));
		StatSig::Init(s_SB->GetDBSize());
		s_SB->Setup();
		s_SB->m_QuerySelf = true;

		vector<uint> TopEvalIdxs;
		P.GetTopEvalIdxs(3, TopEvalIdxs);
		const uint n = SIZE(TopEvalIdxs);
		if (n == 0)
			Die("No evals");
		for (uint k = 0; k < n; ++k)
			{
			ProgressLog("=========================================\n");
			ProgressLog("               FULL %u/%u\n", k+1, n);
			ProgressLog("=========================================\n");

			uint EvalIdx = TopEvalIdxs[k];
			const vector<string> &xv = P.m_xvs[EvalIdx];
			double Area = EvalArea3(xv);
			ProgressLog("FULLDB [%.3g] %s\n", Area, BestVarStr.c_str());
			}
		}

	CloseStdioFile(P.m_fTsv);
	}
