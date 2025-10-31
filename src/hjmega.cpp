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

static void Optimize(
	const string &OptName,
	const vector<string> &SpecLines,
	SCOP40Bench &SB,
	uint LatinBinCount,
	uint HJCount,
	double &Best_y,
	vector<string> &Best_xv)
	{
	Peaker &P = *new Peaker(0, OptName);
	P.Init(SpecLines, EvalArea3);
	s_Peaker = &P;

	s_SB = &SB;

	ProgressLog("=========================================\n");
	ProgressLog("%s latin (%u)\n", OptName.c_str(), LatinBinCount);
	ProgressLog("=========================================\n");
	asserta(LatinBinCount > 0);
	P.RunLatin(LatinBinCount);

	vector<uint> TopEvalIdxs;
	P.GetTopEvalIdxs(3, TopEvalIdxs);
	const uint n = SIZE(TopEvalIdxs);
	if (n == 0)
		Die("No evals %s", OptName.c_str());

	for (uint k = 0; k < n; ++k)
		{
		ProgressLog("=========================================\n");
		ProgressLog("%s HJ %u/%u\n", OptName.c_str(), k+1, n);
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
		ProgressLog("=========================================\n");
		ProgressLog("%s HJ %u/%u converged\n", OptName.c_str(), k+1, n);
		ProgressLog("=========================================\n");
		}
	Best_y = P.m_Best_y;
	Best_xv = P.m_Best_xv;
	ProgressLog("=========================================\n");
	ProgressLog("%s completed\n", OptName.c_str());
	ProgressLog("=========================================\n");
	}

void cmd_hjmega()
	{
	const string SpecFN = g_Arg1;
	asserta(optset_db);
	const string &DBFN = opt(db);

	Peaker::m_fTsv = CreateStdioFile(opt(output2));

	DSSParams::Init(DM_UseCommandLineOption);
	StatSig::Init(SCOP40_DBSIZE);

	SCOP40Bench FullSB;
	FullSB.LoadDB(DBFN);
	FullSB.Setup();

	Log("SpecFN=%s\n", SpecFN.c_str());
	vector<string> SpecLines;
	ReadLinesFromFile(SpecFN, SpecLines);
	for (uint i = 0; i < SIZE(SpecLines); ++i)
		Log("%s\n", SpecLines[i].c_str());

	uint SubsetIters = UINT_MAX;
	uint SubsetPct = UINT_MAX;
	uint LatinCount = UINT_MAX;
	uint HJCount = UINT_MAX;
	string GlobalSpec;
	{
	Peaker TmpP(0, "");
	TmpP.Init(SpecLines, 0);
	GlobalSpec = TmpP.m_GlobalSpec;
	}

	SubsetIters = Peaker::SpecGetInt(GlobalSpec, "sub", UINT_MAX);
	SubsetPct = Peaker::SpecGetInt(GlobalSpec, "subpct", UINT_MAX);
	LatinCount = Peaker::SpecGetInt(GlobalSpec, "latin", UINT_MAX);
	HJCount = Peaker::SpecGetInt(GlobalSpec, "hj", UINT_MAX);
	asserta(SubsetIters != UINT_MAX);
	asserta(LatinCount != UINT_MAX);
	asserta(HJCount != UINT_MAX);

	double Final_y = -1;
	string Final_xss;
	for (uint SubsetIter = 1; SubsetIter <= SubsetIters; ++SubsetIter)
		{
		double Best_y;
		vector<string> Best_xv;
		SCOP40Bench *Subset = new SCOP40Bench;
		FullSB.MakeSubset(*Subset, SubsetPct);
		ProgressLog("Subset %u%%, %u chains\n",
				SubsetPct, Subset->GetDBChainCount());
		string OptName;
		Ps(OptName, "sub%u", SubsetIter);
		Optimize(OptName, SpecLines, *Subset, LatinCount, HJCount,
			Best_y, Best_xv);

		s_SB = &FullSB;	
		string PeakerName;
		Ps(PeakerName, "full%u", SubsetIter);
		Peaker Pfull(0, PeakerName);
		Pfull.Init(SpecLines, EvalArea3);
		Pfull.Evaluate(Best_xv, PeakerName + "_init");
		Pfull.HJ_RunHookeJeeves();
		Pfull.WriteFinalResults(g_fLog);

		if (Pfull.m_Best_y > Final_y)
			{
			Final_y = Pfull.m_Best_y;
			Pfull.xv2xss(Pfull.m_Best_xv, Final_xss);
			}

		delete Subset;
		}

	ProgressLog("\n");
	ProgressLog("FINAL [%.4g] %s\n", Final_y, Final_xss.c_str());

	CloseStdioFile(Peaker::m_fTsv);
	}
