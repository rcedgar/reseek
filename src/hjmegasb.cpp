#include "myutils.h"
#include "statsig.h"
#include "subsetbench.h"
#include "peaker.h"

float SubsetBench_AF_SWFast(
	XDPMem &Mem,
	uint LQ, uint LT,
	float Open, float Ext,
	float * const * SWMx);

static SubsetBench *s_SB;
static Peaker *s_Peaker;

static double EvalSum3(const vector<string> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);
	string VarsStr;
	s_Peaker->xv2xss(xv, VarsStr);
	s_SB->UpdateParamsFromVarStr(VarsStr);
	s_SB->Search(SubsetBench_AF_SWFast, false);
	s_SB->Bench();
	return s_SB->m_Sum3;
	}

static double EvalSum3All(const string &ss)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	s_SB->UpdateParamsFromVarStr(ss);
	s_SB->Search(SubsetBench_AF_SWFast, true);
	s_SB->BenchAll(ss);
	return s_SB->m_Sum3;
	}

static void Optimize(
	const string &OptName,
	const vector<string> &SpecLines,
	SubsetBench &SB,
	double &Best_y,
	string &Best_ss)
	{
	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	uint LatinBinCount = Peaker::SpecGetInt(GlobalSpec, "latin", UINT_MAX);
	uint HJCount = Peaker::SpecGetInt(GlobalSpec, "hj", UINT_MAX);
	asserta(LatinBinCount != UINT_MAX);
	asserta(HJCount != UINT_MAX);

	Peaker &P = *new Peaker(0, OptName);
	P.Init(SpecLines, EvalSum3);
	s_Peaker = &P;
	s_SB = &SB;

	ProgressLog("=========================================\n");
	ProgressLog("%s latin (%u)\n", OptName.c_str(), LatinBinCount);
	ProgressLog("=========================================\n");
	asserta(LatinBinCount > 0);
	P.RunLatin(LatinBinCount);

	vector<uint> TopEvalIdxs;
	P.GetTopEvalIdxs(HJCount, TopEvalIdxs);
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
		string ss;
		Child->xv2xss(Child->m_Best_xv, ss);

		ProgressLog("=========================================\n");
		ProgressLog("%s HJ %u/%u converged: ", OptName.c_str(), k+1, n);
		ProgressLog(" [%.4g] %s\n", Child->m_Best_y, ss.c_str());
		ProgressLog("=========================================\n");

		P.AppendChildResults(*Child);
		delete Child;
		}
	Best_y = P.m_Best_y;
	P.xv2xss(P.m_Best_xv, Best_ss);
	ProgressLog("=========================================\n");
	ProgressLog("%s completed:", OptName.c_str());
	ProgressLog(" [%.4g] %s\n", Best_y, Best_ss.c_str());
	ProgressLog("=========================================\n");
	}

static void Climb(SubsetBench &FullSB, const vector<string> &SpecLines,
	double &Best_y, string &Best_ss)
	{
	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	vector<string> Fields;
	string ParamStr;
	for (uint i = 0; i < SIZE(SpecLines); ++i)
		{
		const string &Line = SpecLines[i];
		if (StartsWith(Line, "#init "))
			{
			ParamStr = Line.substr(6);
			break;
			}
		}
	if (ParamStr.empty())
		Die("Missing #init in spec");

	vector<string> Fields2, VarNames, Init_xv;
	Split(ParamStr, Fields, ';');
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		Split(Fields[i], Fields2, '=');
		asserta(SIZE(Fields2) == 2);
		VarNames.push_back(Fields2[0]);
		Init_xv.push_back(Fields2[1]);
		}
	const uint VarCount = SIZE(VarNames);

	s_SB = &FullSB;	
	string PeakerName;
	Ps(PeakerName, "climb");
	Peaker Pfull(0, PeakerName);
	Pfull.Init(SpecLines, EvalSum3);
	s_Peaker = &Pfull;

	Pfull.Evaluate(Init_xv, PeakerName + "_init");
	Pfull.HJ_RunHookeJeeves();
	Pfull.WriteFinalResults(g_fLog);

	Best_y = Pfull.m_Best_y;
	Pfull.xv2xss(Pfull.m_Best_xv, Best_ss);
	}

void cmd_hjmegasb()
	{
	const string SpecFN = g_Arg1;

	asserta(optset_bspattern);
	asserta(optset_mxpattern);
	asserta(optset_dope);

	const string &DopeFN = opt(dope);
	const string &LookupFN = opt(lookup);
	const string &BSFNPattern = opt(bspattern);
	const string &MxFNPattern = opt(mxpattern);

	Log("SpecFN=%s\n", SpecFN.c_str());
	vector<string> SpecLines;
	ReadLinesFromFile(SpecFN, SpecLines);

	vector<string> VarNames;
	Peaker::GetVarNames(SpecLines, VarNames);
	vector<float> PlaceholderValues(SIZE(VarNames), 1);

	vector<string> AlphaNames;
	vector<string> ScalarNames;
	vector<float> Weights;
	vector<float> ScalarValues;
	SubsetBench::ClassifyParams(
		VarNames, PlaceholderValues, AlphaNames, Weights, ScalarNames, ScalarValues);

	Peaker::m_fTsv = CreateStdioFile(opt(output2));

	SubsetBench SB;
	SB.ReadLookup(LookupFN);
	SB.ReadDope(DopeFN);
	SB.LoadAlphas(AlphaNames, BSFNPattern, MxFNPattern);
	SB.SetWeights(Weights);
	SB.SetScalarParams(ScalarNames, ScalarValues);
	SB.AllocHits();

	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	string Strategy;
	Peaker::SpecGetStr(GlobalSpec, "strategy", Strategy, "");
	if (Strategy == "")
		Die("Missing strategy=");

	double Best_y;
	string Best_ss;
	if (Strategy == "climb")
		Climb(SB, SpecLines, Best_y, Best_ss);
	else if (Strategy == "latinclimb")
		{
		double Best_y;
		vector<string> Best_xv;
		Optimize("LatinClimb", SpecLines, SB, Best_y, Best_ss);
		}
	else
		Die("Bad strategy=%s", Strategy.c_str());
	CloseStdioFile(Peaker::m_fTsv);

	if (optset_topn)
		{
		SB.InitFB();
		vector<uint> Idxs;
		s_Peaker->GetTopEvalIdxs_mindy(opt(topn), 0.01, Idxs);
		const uint n = SIZE(Idxs);
		for (uint i = 0; i < n; ++i)
			{
			uint Idx = Idxs[i];
			string ss;
			const vector<string> &xv = s_Peaker->m_xvs[Idx];
			s_Peaker->xv2xss(xv, ss);
			EvalSum3All(ss);
			}
		}
	}
