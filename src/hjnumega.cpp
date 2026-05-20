#include "myutils.h"
#include "statsig.h"
#include "parabench.h"
#include "peaker.h"
#include "flat_helpers.h"

static ParaBench *s_PB;
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

static void set_params_from_varstr(
	const string &varstr,
	unordered_map<string, float> &name2weight,
	int &IntOpen,
	int &IntExt,
	float &ScaleFactor,
	int &IntSaturatedScore)
	{
	ScaleFactor = 1;
	if (optset_scalef)
		ScaleFactor = float(opt(scalef));
	IntOpen = 0;
	IntExt = 0;
	IntSaturatedScore = 777;
	const uint VarCount = s_Peaker->GetVarCount();
	vector<string> flds;
	Split(varstr, flds, ';');
	asserta(flds.size() == VarCount);
	for (uint VarIdx = 0; VarIdx < VarCount; ++VarIdx)
		{
		const string &name_eq_value = flds[VarIdx];
		vector<string> flds2;
		Split(name_eq_value, flds2, '=');
		asserta(flds2.size() == 2);
		const string &VarName = s_Peaker->GetVarName(VarIdx);
		asserta(flds2[0] == VarName);
		const string &sValue = flds2[1];
		if (VarName == "intopen")
			IntOpen = LocalStrToInt(sValue);
		else if (VarName == "open")
			Die("var=open not supported");
		else if (VarName == "intext")
			IntExt = LocalStrToInt(sValue);
		else if (VarName == "ext")
			Die("var=open not supported");
		else if (VarName == "scale")
			ScaleFactor = (float) StrToFloat(sValue);
		else if (VarName == "gap2")
			Die("var=gap2 not supported");
		else
			name2weight[VarName] = (float) StrToFloat(sValue);
		}
	}

static double EvalSum3(const vector<string> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);

	string varstr;
	s_Peaker->xv2xss(xv, varstr);

	unordered_map<string, float> name2weight;
	int IntOpen;
	int IntExt;
	float Scale;
	int IntSaturatedScore;
	set_params_from_varstr(
		varstr, name2weight, IntOpen, IntExt, Scale, IntSaturatedScore);

	Paralign::set_flat_compound(name2weight,
		Scale, IntOpen, IntExt, IntSaturatedScore);
	s_PB->ClearHitsAndResults();
	s_PB->Search("para", false);
	s_PB->SetScoreOrder();
	s_PB->Bench();
	return s_PB->m_Sum3;
	}

static double EvalSum3_VarStr(ParaBench &PB, const string &VarStr)
	{
	s_PB = &PB;
	unordered_map<string, float> name2weight;
	int IntOpen;
	int IntExt;
	float Scale;
	int IntSaturatedScore;
	set_params_from_varstr(
		VarStr, name2weight, IntOpen, IntExt, Scale, IntSaturatedScore);

	Paralign::set_flat_compound(name2weight,
		Scale, IntOpen, IntExt, IntSaturatedScore);
	s_PB->ClearHitsAndResults();
	s_PB->Search("para", false);
	s_PB->SetScoreOrder();
	s_PB->Bench();
	return s_PB->m_Sum3;
	}

static void Optimize(
	const string &OptName,
	const vector<string> &SpecLines,
	ParaBench &PB,
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
	s_PB = &PB;

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

static void Climb(ParaBench &FullPS, const vector<string> &SpecLines)
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

	s_PB = &FullPS;	
	string PeakerName;
	Ps(PeakerName, "climb");
	Peaker Pfull(0, PeakerName);
	Pfull.Init(SpecLines, EvalSum3);
	s_Peaker = &Pfull;

	Pfull.Evaluate(Init_xv, PeakerName + "_init");
	uint LatinBinCount = Peaker::SpecGetInt(GlobalSpec, "latin", UINT_MAX);
	if (LatinBinCount != UINT_MAX)
		Pfull.RunLatin(LatinBinCount);
	Pfull.HJ_RunHookeJeeves();
	Pfull.WriteFinalResults(g_fLog);
	}

static void SubClimb(ParaBench &FullPS, const vector<string> &SpecLines)
	{
	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	uint SubsetIters = Peaker::SpecGetInt(GlobalSpec, "sub", UINT_MAX);
	uint SubsetPct = Peaker::SpecGetInt(GlobalSpec, "subpct", UINT_MAX);
	asserta(SubsetIters != UINT_MAX);

	double Final_y = -1;
	string Final_xss;

	ParaBench &Subset = *new ParaBench;
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

		s_PB = &FullPS;	
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
	asserta(!optset_scale);

	const string SpecFN = g_Arg1;
	Log("SpecFN=%s\n", SpecFN.c_str());
	vector<string> SpecLines;
	ReadLinesFromFile(SpecFN, SpecLines);

	vector<string> alpha_names;
	get_alpha_names_from_peaker_spec_file_lines(
		SpecLines, alpha_names);

	const uint AlphaCount = SIZE(alpha_names);
	asserta(AlphaCount > 0);
	vector<float> weights(AlphaCount, 1.0f); // placeholder

	flat_params::init_from_alphadir(opt(alphadir), alpha_names);
	asserta(flat_params::m_nfeat == AlphaCount);

	asserta(optset_db);
	const string &DBFN = opt(db);

	//void OpenOutputFiles();
	//OpenOutputFiles();
	Peaker::m_fTsv = CreateStdioFile(opt(output2));

	ParaBench FullPS;
	FullPS.GetByteSeqs(DBFN, "nuletters");
	FullPS.SetLookupFromLabels();
	FullPS.to_hexfasta(opt(hexfasta));
	if (optset_varstr)
		{
		Peaker &P = *new Peaker(0, "varstr");
		P.Init(SpecLines, EvalSum3);
		s_Peaker = &P;
		EvalSum3_VarStr(FullPS, opt(varstr));
		FullPS.WriteHits(opt(output));
		Paralign::LogMatrix();
		return;
		}

	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	string Strategy;
	Peaker::SpecGetStr(GlobalSpec, "strategy", Strategy, "");
	if (Strategy == "")
		Die("Missing strategy=");

	if (Strategy == "climb")
		{
		Climb(FullPS, SpecLines);
		}
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
	Paralign::LogMatrix();
	}
