#include "myutils.h"
#include "statsig.h"
#include "parasearch.h"
#include "peaker.h"

static ParaSearch *s_PS;
static Peaker *s_Peaker;

static void GetFeaturesFromVarNames(const Peaker &P, vector<FEATURE> &Fs)
	{
	const uint VarCount = P.GetVarCount();
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		FEATURE F = StrToFeature(P.m_VarNames[VarIdx].c_str(), true);
		if (F != FEATURE(UINT_MAX))
			Fs.push_back(F);
		}
	}

static int LocalStrToInt(const string &s)
	{
	float f = StrToFloatf(s);
	int i = int(round(f));
	float f2 = float(i);
	asserta(f == f2);
	return i;
	}

static double EvalSum3(const vector<string> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);
	string VarsStr;

	int ScaleFactor = 1;
	if (optset_scale)
		ScaleFactor = opt(scale);
	int Open = 0;
	int Ext = 0;
	int SaturatedScore = 777;
	vector<float> Weights;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		string sValue = xv[VarIdx];
		const string &VarName = s_Peaker->GetVarName(VarIdx);
		if (VarName == "open")
			Open = LocalStrToInt(sValue);
		else if (VarName == "ext")
			Ext = LocalStrToInt(sValue);
		else if (VarName == "scale")
			Die("ver=scale not supported");
		else
			{
			FEATURE F = StrToFeature(VarName.c_str());
			asserta(F == ParaSearch::m_NuFs[SIZE(Weights)]);
			Weights.push_back(StrToFloatf(sValue));
			}
		}

	Paralign::SetCompoundMx(ParaSearch::m_NuFs, Weights,
		ScaleFactor, Open, Ext, SaturatedScore);
	s_PS->ClearHitsAndResults();
	s_PS->Search("para", false);
	s_PS->SetScoreOrder();
	s_PS->Bench();
	return s_PS->m_Sum3;
	}

static double EvalSum3_VarStr(ParaSearch &PS, const string &VarStr)
	{
	vector<string> Fields, Fields2;
	Split(VarStr, Fields, ';');
	const uint VarCount = SIZE(Fields);

	int ScaleFactor = 1;
	if (optset_scale)
		ScaleFactor = opt(scale);
	int Open = 0;
	int Ext = 0;
	int SaturatedScore = 777;
	vector<float> Weights;
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const string &Field = Fields[VarIdx];
		Split(Field, Fields2, '=');
		asserta(SIZE(Fields2) == 2);
		const string &VarName = Fields2[0];
		const string &sValue = Fields2[1];
		if (VarName == "open")
			Open = StrToInt(sValue);
		else if (VarName == "ext")
			Ext = StrToInt(sValue);
		else
			{
			FEATURE F = StrToFeature(VarName.c_str());
			asserta(F == ParaSearch::m_NuFs[SIZE(Weights)]);
			Weights.push_back(StrToFloatf(sValue));
			}
		}

	Paralign::SetCompoundMx(ParaSearch::m_NuFs,
		Weights, ScaleFactor, Open, Ext, SaturatedScore);
	PS.ClearHitsAndResults();
	PS.Search("para", false);
	PS.SetScoreOrder();
	PS.Bench();
	return PS.m_Sum3;
	}

static void Optimize(
	const string &OptName,
	const vector<string> &SpecLines,
	ParaSearch &PS,
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
	P.Init(SpecLines, EvalSum3);
	s_Peaker = &P;
	s_PS = &PS;

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

static void Climb(ParaSearch &FullPS, const vector<string> &SpecLines)
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

	s_PS = &FullPS;	
	string PeakerName;
	Ps(PeakerName, "climb");
	Peaker Pfull(0, PeakerName);
	Pfull.Init(SpecLines, EvalSum3);
	//asserta(Pfull.GetVarCount() == VarCount);
	//asserta(Pfull.m_VarNames == VarNames);
	s_Peaker = &Pfull;

	Pfull.Evaluate(Init_xv, PeakerName + "_init");
	Pfull.HJ_RunHookeJeeves();
	Pfull.WriteFinalResults(g_fLog);
	}

static void SubClimb(ParaSearch &FullPS, const vector<string> &SpecLines)
	{
	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	uint SubsetIters = Peaker::SpecGetInt(GlobalSpec, "sub", UINT_MAX);
	uint SubsetPct = Peaker::SpecGetInt(GlobalSpec, "subpct", UINT_MAX);
	asserta(SubsetIters != UINT_MAX);

	double Final_y = -1;
	string Final_xss;

	ParaSearch &Subset = *new ParaSearch;
	for (uint SubsetIter = 1; SubsetIter <= SubsetIters; ++SubsetIter)
		{
		double Best_y;
		vector<string> Best_xv;
		FullPS.MakeSubset(Subset, SubsetPct);
		ProgressLog("Subset %u%%, %u chains\n",
			SubsetPct, SIZE(Subset.m_Labels));
		string OptName;
		Ps(OptName, "sub%u", SubsetIter);
		Optimize(OptName, SpecLines, Subset, Best_y, Best_xv);

		s_PS = &FullPS;	
		string PeakerName;
		Ps(PeakerName, "all%u", SubsetIter);
		Peaker Pfull(0, PeakerName);
		Pfull.Init(SpecLines, EvalSum3);
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
	ProgressLog("FINAL subclimb [%.4g] %s\n",
		Final_y, Final_xss.c_str());
	Log("@TSV@");
	Log("\t%.4g", Final_y);
	Log("\t%s", Final_xss.c_str());
	Log("\n");
	}

void cmd_hjnumega()
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
	ParaSearch::m_NuFs = Fs;

	StatSig::Init(SCOP40_DBSIZE);

	asserta(optset_db);
	const string &DBFN = opt(db);

	void OpenOutputFiles();
	OpenOutputFiles();
	Peaker::m_fTsv = CreateStdioFile(opt(output2));

	//DSSParams::Init(DM_UseCommandLineOption);
	//DSSParams::OverwriteFeatures(Fs, Weights);

	ParaSearch FullPS;
	FullPS.GetByteSeqs(DBFN, "nuletters");
	FullPS.SetLookupFromLabels();
	if (optset_varstr)
		{
		EvalSum3_VarStr(FullPS, opt(varstr));
		return;
		}

	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	string Strategy;
	Peaker::SpecGetStr(GlobalSpec, "strategy", Strategy, "");
	if (Strategy == "")
		Die("Missing strategy=");

	if (Strategy == "climb")
		Climb(FullPS, SpecLines);
	else if (Strategy == "subclimb")
		SubClimb(FullPS, SpecLines);
	else if (Strategy == "latinclimb")
		{
		double Best_y;
		vector<string> Best_xv;
		Optimize("LatinClimb", SpecLines, FullPS, Best_y, Best_xv);
		}
	else
		Die("Bad strategy=%s", Strategy.c_str());

	CloseStdioFile(Peaker::m_fTsv);
	}
