#include "myutils.h"
#include "flat_bench.h"
#include "peaker.h"

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values);

static flat_bench *s_FB;
static Peaker *s_Peaker;

static double EvalSum3(const vector<string> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);
	string VarStr;
	s_Peaker->xv2xss(xv, VarStr);
	s_FB->UpdateParamsFromVarStr(VarStr);
	s_FB->ClearHitsAndResults();
	uint ThreadCount = GetRequestedThreadCount();
	s_FB->Search(ThreadCount, false, optset_dope, UINT_MAX);
	s_FB->SetScoreOrder_Parallel();
	s_FB->Bench();
	return s_FB->m_Sum3;
	}

static void Optimize(
	const string &OptName,
	const vector<string> &SpecLines,
	flat_bench &FB,
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
	s_FB = &FB;

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
	ProgressLog("FINAL %.3g\n", Best_y);
	ProgressLog("=========================================\n");
	}

static void Climb(flat_bench &FullFB, const vector<string> &SpecLines)
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

	s_FB = &FullFB;	
	string PeakerName;
	Ps(PeakerName, "climb");
	Peaker Pfull(0, PeakerName);
	Pfull.Init(SpecLines, EvalSum3);
	s_Peaker = &Pfull;

	Pfull.Evaluate(Init_xv, PeakerName + "_init");
	Pfull.HJ_RunHookeJeeves();
	Pfull.WriteFinalResults(g_fLog);
	}

static void Resume(flat_bench &FullFB, const vector<string> &SpecLines)
	{
	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	s_FB = &FullFB;	
	string PeakerName;
	Ps(PeakerName, "climb");
	Peaker Pfull(0, PeakerName);
	Pfull.Init(SpecLines, EvalSum3);
	Pfull.LoadTSV(opt(input2));
	s_Peaker = &Pfull;

	Pfull.HJ_RunHookeJeeves();
	Pfull.WriteFinalResults(g_fLog);
	}

static void SubClimb(
	flat_bench &FullFB, 
	flat_bench &SubsetFB, 
	const vector<string> &SpecLines)
	{
	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	uint SubsetIters = Peaker::SpecGetInt(GlobalSpec, "sub", UINT_MAX);
	asserta(SubsetIters != UINT_MAX);

	double Final_y = -1;
	string Final_xss;

	for (uint SubsetIter = 1; SubsetIter <= SubsetIters; ++SubsetIter)
		{
		double Best_y;
		vector<string> Best_xv;
		ProgressLog("Subset %u chains\n", SubsetFB.m_fp.get_nprof());
		string OptName;
		Ps(OptName, "sub%u", SubsetIter);
		Optimize(OptName, SpecLines, SubsetFB, Best_y, Best_xv);

		s_FB = &FullFB;	
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
	ProgressLog("FINAL subclimb [%.4g] %s\n", Final_y, Final_xss.c_str());
	Log("@TSV@");
	Log("\t%.4g", Final_y);
	Log("\t%s", Final_xss.c_str());
	Log("\n");
	}

void get_feature_names_from_peaker_spec_file_lines(
	vector<string> &lines,
	vector<string> &alpha_names)
	{
	alpha_names.clear();
	vector<string> flds;
	for (size_t i = 0; i < lines.size(); ++i)
		{
		const string &line = lines[i];
		if (!StartsWith(line, "var="))
			continue;
		Split(line, flds, ';');
		const string var_eq_name = flds[0];
		Split(var_eq_name, flds, '=');
		asserta(flds.size() == 2);
		const string &name = flds[1];
		if (line.find("isalpha=yes;") != string::npos)
			alpha_names.push_back(name);
		}
	}

void cmd_flat_hjmega()
	{
	asserta(optset_fapattern);
	asserta(optset_mxpattern);
	asserta(!optset_varstr); // must assert agrees with spec

	const string SpecFN = g_Arg1;
	Log("SpecFN=%s\n", SpecFN.c_str());
	vector<string> SpecLines;
	ReadLinesFromFile(SpecFN, SpecLines);

	vector<string> AlphaNames;
	get_feature_names_from_peaker_spec_file_lines(
		SpecLines, AlphaNames);

	void OpenOutputFiles();
	OpenOutputFiles();
	Peaker::m_fTsv = CreateStdioFile(opt(output2));

	flat_bench FullFB;
	FullFB.ReadLookup(opt(lookup));
	flat_features::load_alphas(AlphaNames, opt(mxpattern));
	FullFB.load_profiles(opt(fapattern));
	vector<flat_chain_t *> chains;
	read_flat_chains(opt(input), chains);
	FullFB.set_distmxs(chains);
	if (optset_dope)
		FullFB.ReadDope(opt(dope));
	FullFB.Alloc();

	if (optset_input2)
		{
		Resume(FullFB, SpecLines);
		CloseStdioFile(Peaker::m_fTsv);
		return;
		}

	string GlobalSpec;
	Peaker::GetGlobalSpec(SpecLines, GlobalSpec);

	string Strategy;
	Peaker::SpecGetStr(GlobalSpec, "strategy", Strategy, "");

	if (Strategy == "")
		Die("Missing strategy=");

	if (Strategy == "climb")
		Climb(FullFB, SpecLines);
	else if (Strategy == "subclimb")
		{
		asserta(optset_subdope);
		asserta(optset_sublookup);

		flat_bench SubsetFB;
		SubsetFB.ReadLookup(opt(sublookup));
		SubsetFB.load_profiles(opt(fapattern));
		SubsetFB.set_distmxs(chains);
		SubsetFB.LogParams();
		SubsetFB.ReadDope(opt(subdope));
		SubsetFB.Alloc();
		SubClimb(FullFB, SubsetFB, SpecLines);
		}
	else if (Strategy == "latinclimb")
		{
		double Best_y;
		vector<string> Best_xv;
		Optimize("LatinClimb", SpecLines, FullFB, Best_y, Best_xv);
		}
	else
		Die("Bad strategy=%s", Strategy.c_str());

	CloseStdioFile(Peaker::m_fTsv);
	}
