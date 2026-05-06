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

void cmd_flat_hjmega_bracket()
	{
	asserta(optset_fapattern);
	asserta(optset_mxpattern);

	const string &VarStr = g_Arg1;

	vector<string> param_names;
	vector<float> param_values;
	ParseVarStr(VarStr, param_names, param_values);

	vector<string> alpha_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		alpha_names, weights,
		scalar_names, scalar_values);

	void OpenOutputFiles();
	OpenOutputFiles();
	Peaker::m_fTsv = CreateStdioFile(opt(output2));

	flat_bench FullFB;
	FullFB.ReadLookup(opt(lookup));
	flat_features::load_alphas(alpha_names, opt(mxpattern));
	FullFB.load_profiles(opt(fapattern));
	vector<flat_chain_t *> chains;
	read_flat_chains(opt(input), chains);
	FullFB.set_distmxs(chains);
	if (optset_dope)
		FullFB.ReadDope(opt(dope));
	FullFB.Alloc();

	vector<string> spec_lines;
	const uint nalpha = uint(alpha_names.size());
	const uint nscalar = uint(scalar_names.size());
	const uint nvar = nalpha + nscalar;
	for (uint i = 0; i < nscalar; ++i)
		{
		string line;
		// var=dali;min=0.00655;max=0.00792;sigfig=2;
		Ps(line, "var=%s;weight=no;alpha=no;sigfig=5;", scalar_names[i].c_str());
		spec_lines.push_back(line);
		}
	for (uint i = 0; i < nalpha; ++i)
		{
		string line;
		// var=aa20;min=0.373;max=0.451;sigfig=2;weight=yes;isalpha=yes;
		Ps(line, "var=%s;weight=yes;alpha=yes;sigfig=5;", alpha_names[i].c_str());
		spec_lines.push_back(line);
		}

	s_FB = &FullFB;	
	string PeakerName;
	Ps(PeakerName, "bracket");
	Peaker Pfull(0, PeakerName);
	Pfull.Init(spec_lines, EvalSum3);
	s_Peaker = &Pfull;
	vector<string> flds;

	Pfull.Evaluate(VarStr, "init");

	vector<string> plus_value_strs;
	vector<string> minus_value_strs;
	vector<double> plus_ys;
	vector<double> minus_ys;
	vector<double> rates(nvar, FLT_MAX);

	for (;;)
		{
		double saved_besty = Pfull.m_Best_y;
		Pfull.Bracket_FindNeighbors(
			0.02, 0.2, 5, 
			plus_value_strs, minus_value_strs,
			plus_ys, minus_ys,
			rates);
		asserta(plus_value_strs.size() == nvar);
		double y = Pfull.m_Best_y;
		if (y == saved_besty)
			{
			ProgressLog("\nCONVERGED\n");
			break;
			}
		double pct = GetPct(y - saved_besty, saved_besty);
		ProgressLog("\nIMPROVED y %.5g => %.5g (%.2f%%)\n",
			saved_besty, y, pct);
		}

	double besty = Pfull.m_Best_y;

	ProgressLog("\n\n");
	ProgressLog("FINAL BRACKET\n\n");
	ProgressLog("%12.12s", "var");
	ProgressLog("  %12.12s", "+value");
	ProgressLog("  %9.9s", "+dy%");
	ProgressLog("  %12.12s", "-value");
	ProgressLog("  %9.9s", "-dy%");
	ProgressLog("  %12.12s", "rate");
	ProgressLog("\n");
	for (uint VarIdx = 0; VarIdx < nvar; ++VarIdx)
		{
		double plusy = plus_ys[VarIdx];
		double minusy = minus_ys[VarIdx];
		double plusy_pct = GetPct(besty - plusy, besty);
		double minusy_pct = GetPct(besty - minusy, besty);
		ProgressLog("%12.12s", Pfull.GetVarName(VarIdx));
		ProgressLog("  %12.12s", plus_value_strs[VarIdx]);
		ProgressLog("  %8.2f%%", plusy_pct);
		ProgressLog("  %12.12s", minus_value_strs[VarIdx]);
		ProgressLog("  %8.2f%%", minusy_pct);
		ProgressLog("  %10.5g", rates[VarIdx]);
		ProgressLog("\n");
		}
	Pfull.WriteFinalResults(g_fLog);

	CloseStdioFile(Peaker::m_fTsv);
	}
