#include "myutils.h"
#include "statsig.h"
#include "scop40bench.h"
#include "peaker.h"

/***
(+) Random exploration like Latin hypercube quickly reaches diminishing returns
(+) Careful hill-climbing from a good start is essential
(+) Preferred strategy per 2025-11-05 is "subclimb" (X is vector of parameter values):
		Repeat 3 times (i):
			Select random subset 25% of SCOP40[c] chains ("subpct=25;" in global sec),
				Eval of one X 16x faster than working with full set
			Perform Latin hypercube on subset
			Take top 3 points ("hj=3;" in global spec) and hill-climb (HJ)
			From these 3, save point X_i with best y
			Hill-climb (HJ) full set starting at X_i
		Final X is best of the 3 full HJs
(*) Subclimb time with ~10 params ~20 mins on a/b subset, 2 - 4 hours on SCOP40[c]
(*) Don't need open+ext, gap2 works fine (ext = open/10)
(*) -raw avoids calculating self-rev scores & tuning E-value parameters,
	  gives similar parameters so this is the way to choose alphabet components.
(*) The three lowest-weight alphabets in v2.7 consistently optimize to weight=0
		AddFeature(FEATURE_DstNxtHlx,	0.00475462f);
		AddFeature(FEATURE_StrandDens,	0.0183853f);
		AddFeature(FEATURE_NormDens,	0.00384384f);
***/

static SCOP40Bench *s_SB;
static Peaker *s_Peaker;

static void AssertProfileSize(const DBSearcher &DBS, uint FeatureCount)
	{
	const uint ChainCount = DBS.GetDBChainCount();
	for (uint i = 0; i < ChainCount; ++i)
		{
		const vector<vector<byte> > &Profile = *DBS.m_DBProfiles[i];
		asserta(SIZE(Profile) == FeatureCount);
		}
	}

static void GetFeaturesFromVarNames(const Peaker &P, vector<FEATURE> &Fs)
	{
	DSSParams::m_Features.clear();
	DSSParams::m_Weights.clear();

	const uint VarCount = P.GetVarCount();
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		FEATURE F = StrToFeature(P.m_VarNames[VarIdx].c_str(), true);
		if (F != FEATURE(UINT_MAX))
			Fs.push_back(F);
		}
	}

static double EvalArea0(const vector<string> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);
	string VarsStr;
	s_Peaker->xv2xss(xv, VarsStr);
	DSSParams::SetParamsFromStr(VarsStr);
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
	DSSParams::SetParamsFromStr(VarsStr);
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
	string ParamsStr;
	DSSParams::GetParamsStr(ParamsStr);
	ProgressLog("ParamsStr:%s\n", ParamsStr.c_str());

	void OpenOutputFiles();
	OpenOutputFiles();

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
	double &Best_y,
	vector<string> &Best_xv)
	{
	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	uint LatinBinCount = Peaker::SpecGetInt(GlobalSpec, "latin", UINT_MAX);
	uint HJCount = Peaker::SpecGetInt(GlobalSpec, "hj", UINT_MAX);
	asserta(LatinBinCount != UINT_MAX);
	asserta(HJCount != UINT_MAX);

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

static void Climb(SCOP40Bench &FullSB, const vector<string> &SpecLines)
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
	Pfull.Init(SpecLines, EvalArea3);
	//asserta(Pfull.GetVarCount() == VarCount);
	//asserta(Pfull.m_VarNames == VarNames);
	s_Peaker = &Pfull;

	Pfull.Evaluate(Init_xv, PeakerName + "_init");
	Pfull.HJ_RunHookeJeeves();
	Pfull.WriteFinalResults(g_fLog);
	}

static void SubClimb(SCOP40Bench &FullSB, const vector<string> &SpecLines)
	{
	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	uint SubsetIters = Peaker::SpecGetInt(GlobalSpec, "sub", UINT_MAX);
	uint SubsetPct = Peaker::SpecGetInt(GlobalSpec, "subpct", UINT_MAX);
	asserta(SubsetIters != UINT_MAX);

	double Final_y = -1;
	string Final_xss;

	// DBSearcher d'tor crashes with subset
	SCOP40Bench &Subset = *new SCOP40Bench;
	for (uint SubsetIter = 1; SubsetIter <= SubsetIters; ++SubsetIter)
		{
		double Best_y;
		vector<string> Best_xv;
		FullSB.MakeSubset(Subset, SubsetPct);
		ProgressLog("Subset %u%%, %u chains\n",
				SubsetPct, Subset.GetDBChainCount());
		string OptName;
		Ps(OptName, "sub%u", SubsetIter);
		Optimize(OptName, SpecLines, Subset, Best_y, Best_xv);

		s_SB = &FullSB;	
		string PeakerName;
		Ps(PeakerName, "all%u", SubsetIter);
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
		}

	ProgressLog("\n");
	ProgressLog("FINAL subclimb [%.4g] %s\n", Final_y, Final_xss.c_str());
	}

void cmd_hjmega()
	{
	const string SpecFN = g_Arg1;
	Log("SpecFN=%s\n", SpecFN.c_str());
	vector<string> SpecLines;
	ReadLinesFromFile(SpecFN, SpecLines);

	// Just to get features
	Peaker Ptmp(0, "tmp");
	Ptmp.Init(SpecLines, 0);
	vector<FEATURE> Fs;
	GetFeaturesFromVarNames(Ptmp, Fs);
	const uint FeatureCount = SIZE(Fs);
	asserta(FeatureCount > 0);
	vector<float> Weights;
	for (uint i = 0; i < FeatureCount; ++i)
		Weights.push_back(1);

	StatSig::Init(SCOP40_DBSIZE);

	asserta(optset_db);
	const string &DBFN = opt(db);

	void OpenOutputFiles();
	OpenOutputFiles();
	Peaker::m_fTsv = CreateStdioFile(opt(output2));

	DSSParams::Init(DM_UseCommandLineOption);
	DSSParams::OverwriteFeatures(Fs, Weights);

	SCOP40Bench FullSB;
	FullSB.LoadDB(DBFN);
	FullSB.Setup();
	FullSB.m_RecalcSelfRevScores = true;
	AssertProfileSize(FullSB, FeatureCount);

	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	string Strategy;
	Peaker::SpecGetStr(GlobalSpec, "strategy", Strategy, "");
	if (Strategy == "")
		Die("Missing strategy=");

	if (Strategy == "climb")
		Climb(FullSB, SpecLines);
	else if (Strategy == "subclimb")
		SubClimb(FullSB, SpecLines);
	else if (Strategy == "latinclimb")
		{
		double Best_y;
		vector<string> Best_xv;
		Optimize("LatinClimb", SpecLines, FullSB, Best_y, Best_xv);
		}
	else
		Die("Bad strategy=%s", Strategy.c_str());

	CloseStdioFile(Peaker::m_fTsv);
	}
