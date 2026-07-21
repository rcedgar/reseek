#include "myutils.h"
#include "flat_bench2.h"
#include "flat_helpers.h"
#include "peaker.h"
#include "paralign.h"
#include "chain_data.h"
#include <unordered_map>

static flat_bench2 *s_FB;
static Peaker *s_Peaker;
static flat_params *s_params = 0;

static double EvalSum3(const vector<string> &xv)
	{
	asserta(s_Peaker != 0);
	const uint VarCount = s_Peaker->GetVarCount();
	asserta(SIZE(xv) == VarCount);
	string VarStr;
	s_Peaker->xv2xss(xv, VarStr);
	vector<string> param_names;
	vector<float> param_values;
	parse_varstr(VarStr, param_names, param_values);
	s_FB->update_params(param_names, param_values);
	s_FB->ClearHitsAndResults();
	uint ThreadCount = GetRequestedThreadCount();
	s_FB->search(ThreadCount, false);
	s_FB->SetScoreOrder_Parallel();
	double Sum3 = s_FB->Bench();
	if (opt(sffp)) return s_FB->m_SFFP;
	return Sum3;
	}

static void Optimize(
	const string &OptName,
	const vector<string> &SpecLines,
	flat_bench2 &FB,
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

static void Climb(flat_bench2 &FullFB, const vector<string> &SpecLines)
	{
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

	vector<string> Fields, Fields2;
	unordered_map<string, string> NameToInit;
	Split(ParamStr, Fields, ';');
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		if (Fields[i].empty())
			continue;
		Split(Fields[i], Fields2, '=');
		asserta(SIZE(Fields2) == 2);
		NameToInit[Fields2[0]] = Fields2[1];
		}

	s_FB = &FullFB;
	string PeakerName;
	Ps(PeakerName, "climb");
	Peaker Pfull(0, PeakerName);
	Pfull.Init(SpecLines, EvalSum3);
	s_Peaker = &Pfull;

	const uint VarCount = Pfull.GetVarCount();
	vector<string> Init_xv;
	for (uint i = 0; i < VarCount; ++i)
		{
		const string &name = Pfull.m_VarNames[i];
		unordered_map<string, string>::const_iterator it = NameToInit.find(name);
		if (it == NameToInit.end())
			Die("Climb #init missing peaker var '%s'", name.c_str());
		Init_xv.push_back(it->second);
		}

	Pfull.Evaluate(Init_xv, PeakerName + "_init");
	Pfull.HJ_RunHookeJeeves();
	Pfull.WriteFinalResults(g_fLog);
	}

static void Resume(flat_bench2 &FullFB, const vector<string> &SpecLines)
	{
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

static void DoConst(flat_bench2 &FullFB, const vector<string> &SpecLines)
	{
	s_FB = &FullFB;
	string PeakerName;
	Ps(PeakerName, "const");
	Peaker Pfull(0, PeakerName);
	Pfull.Init(SpecLines, EvalSum3);
	s_Peaker = &Pfull;
	vector<string> xv;
	Pfull.GetAllConst_xv(xv);
	Pfull.Evaluate(xv, PeakerName + "_const");
	Pfull.WriteFinalResults(g_fLog);
	}

static void SubClimb(
	flat_bench2 &FullFB,
	flat_bench2 &SubsetFB,
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
		ProgressLog("Subset %u domains\n", SubsetFB.m_look->get_ndom());
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

void get_alpha_names_from_peaker_spec_file_lines(
	vector<string> &lines,
	vector<string> &alpha_names)
	{
	validate_flat_peaker_spec(lines, alpha_names);
	}

static void setup_bench(
	flat_bench2 &FB,
	flat_params &params,
	const vector<flat_chain_t *> &chains,
	const string &dope_fn,
	bool have_dope)
	{
	FB.m_params = &params;
	FB.load_chains(chains);
	chain_data::log_mem_stats(params, FB.m_cdvec, FB.m_look->get_ndom());
	if (have_dope)
		FB.ReadDope(dope_fn);
	FB.Alloc();
	}

void cmd_flat_hjmega()
	{
	asserta(!optset_fapattern);
	asserta(!optset_mxpattern);
	asserta(!optset_varstr); // must assert agrees with spec

	const string SpecFN = g_Arg1;
	Log("SpecFN=%s\n", SpecFN.c_str());
	vector<string> SpecLines;
	ReadLinesFromFile(SpecFN, SpecLines);

	vector<string> alpha_names;
	validate_flat_peaker_spec(SpecLines, alpha_names);

	void OpenOutputFiles();
	OpenOutputFiles();
	Peaker::m_fTsv = CreateStdioFile(opt(output2));

	Paralign::set_final_nu();

	flat_params params;
	const string alphadir = optset_alphadir ? string(opt(alphadir)) : "";
	params.init_from_alphadir(alphadir, alpha_names);
	s_params = &params;

	vector<flat_chain_t *> chains;
	read_flat_chains(opt(input), chains);

	flat_bench2 FullFB;
	FullFB.ReadLookup(opt(lookup));
	setup_bench(FullFB, params, chains, opt(dope), optset_dope);

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

	if (Strategy == "const")
		DoConst(FullFB, SpecLines);
	else if (Strategy == "climb")
		Climb(FullFB, SpecLines);
	else if (Strategy == "subclimb")
		{
		asserta(optset_subdope);
		asserta(optset_sublookup);

		flat_bench2 SubsetFB;
		SubsetFB.ReadLookup(opt(sublookup));
		setup_bench(SubsetFB, params, chains, opt(subdope), true);
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
