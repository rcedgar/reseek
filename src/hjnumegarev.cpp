#include "myutils.h"
#include "statsig.h"
#include "parasearch.h"
#include "peaker.h"

static ParaSearch *s_PS;
static Peaker *s_Peaker;

static void GetFeatures(const string &s,
	vector<FEATURE> &Fs, vector<float> &Weights)
	{
	Fs.clear();
	Weights.clear();
	vector<string> Fields;
	Split(s, Fields, ';');
	const uint n = SIZE(Fields);
	for (uint Idx = 0; Idx < n; ++Idx)
		{
		vector<string> Fields2;
		Split(Fields[Idx], Fields2, '=');
		asserta(SIZE(Fields2) == 2);
		FEATURE F = StrToFeature(Fields2[0].c_str());
		float Weight = StrToFloatf(Fields2[1]);
		Fs.push_back(F);
		Weights.push_back(Weight);
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
	asserta(VarCount == 2);
	asserta(s_Peaker->m_VarNames[0] == "selfw");
	asserta(s_Peaker->m_VarNames[1] == "revw");

	float SelfWeight = StrToFloatf(xv[0]);
	float RevWeight = StrToFloatf(xv[1]);
	s_PS->BenchRev("EvalSum3()", SelfWeight, RevWeight);
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
	PS.Search("para", true);
	PS.SetScoreOrder();
	PS.Bench();
	return PS.m_Sum3;
	}

static void Optimize(
	const vector<string> &SpecLines,
	ParaSearch &PS,
	double &Best_y,
	vector<string> &Best_xv)
	{
	const string OptName("latinclimb");
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

static void Climb(ParaSearch &PS, const vector<string> &SpecLines)
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

	s_PS = &PS;	
	string PeakerName;
	Ps(PeakerName, "climb");
	Peaker Pfull(0, PeakerName);
	Pfull.Init(SpecLines, EvalSum3);
	s_Peaker = &Pfull;

	Pfull.Evaluate(Init_xv, PeakerName + "_init");
	Pfull.HJ_RunHookeJeeves();
	Pfull.WriteFinalResults(g_fLog);
	}

void cmd_hjnumegarev()
	{
	vector<FEATURE> Fs;
	vector<float> Weights;

	GetFeatures(g_Arg1, Fs, Weights);
	const uint FeatureCount = SIZE(Fs);
	asserta(FeatureCount > 0);

	asserta(!optset_alignmethod);
	asserta(optset_db);

	asserta(optset_db);
	const string &DBFN = opt(db);

	asserta(optset_intopen);
	asserta(optset_intext);
	asserta(optset_scale);
	int Scale = opt(scale);
	int IntOpen = opt(intopen);
	int IntExt = opt(intext);
	Paralign::SetCompoundMx(Fs, Weights, Scale, 
		IntOpen, IntExt, 777);

	ParaSearch PS;
	ParaSearch::m_NuFs = Fs;
	PS.GetByteSeqs(DBFN, "nuletters");
	PS.SetLookupFromLabels();
	PS.m_DoReverse = true;
	PS.SetSelfScores_rev("para");
	PS.Search("para", true);
	PS.SetScoreOrder();
	PS.WriteRevTsv(opt(output2));
	PS.Bench("Bench()");
	PS.BenchRev("BenchRev(0, 0)", 0, 0);

	vector<string> SpecLines;
	SpecLines.push_back("latin=32;");
	SpecLines.push_back("rates=1.3,1.05,1.02;");
	SpecLines.push_back("hj=2;");
	SpecLines.push_back("var=selfw;min=0;max=1;");
	SpecLines.push_back("var=revw;min=0;max=1;");

	double Best_y;
	vector<string> Best_xv;
	Optimize(SpecLines, PS, Best_y, Best_xv);
	}
