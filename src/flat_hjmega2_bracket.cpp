#include "myutils.h"
#include "flat_bench2.h"
#include "flat_helpers.h"
#include "peaker.h"

static flat_bench2 *s_FB;
static Peaker *s_Peaker;

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
	return Sum3;
	}

void cmd_flat_hjmega2_bracket()
	{
	asserta(optset_varstr);
	const string &varstr = opt(varstr);
	vector<string> param_names;
	vector<float> param_values;
	parse_varstr(opt(varstr), param_names, param_values);

	vector<string> alpha_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_classify_params(
		param_names, param_values,
		alpha_names, weights,
		scalar_names, scalar_values);

	const string &alphadir = opt(alphadir);
	flat_params params;
	params.init_from_alphadir(alphadir, alpha_names);
	params.logme();

	vector<flat_chain_t *> chains;
	read_flat_chains(g_Arg1, chains);

	flat_bench2 FB;
	FB.m_params = &params;
	FB.ReadLookup(opt(lookup));
	FB.load_chains(chains);
	FB.update_params(param_names, param_values);


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

	s_FB = &FB;	
	string PeakerName;
	Ps(PeakerName, "bracket");
	Peaker Pfull(0, PeakerName);
	Pfull.Init(spec_lines, EvalSum3);
	s_Peaker = &Pfull;
	vector<string> flds;

	Pfull.Evaluate(varstr, "init");

	vector<string> plus_value_strs;
	vector<string> minus_value_strs;
	vector<double> plus_ys;
	vector<double> minus_ys;
	vector<double> rates(nvar, FLT_MAX);

	const double mindypct = (optset_mindypct ? opt(mindypct) : 0.02f);
	const double maxdypct = (optset_mindypct ? opt(maxdypct) : 0.2f);

	for (;;)
		{
		double saved_besty = Pfull.m_Best_y;
		Pfull.Bracket_FindNeighbors(
			mindypct, maxdypct, 5, 
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
		ProgressLog("  %12.12s", plus_value_strs[VarIdx].c_str());
		ProgressLog("  %8.2f%%", plusy_pct);
		ProgressLog("  %12.12s", minus_value_strs[VarIdx].c_str());
		ProgressLog("  %8.2f%%", minusy_pct);
		ProgressLog("  %10.5g", rates[VarIdx]);
		ProgressLog("\n");
		}
	Pfull.WriteFinalResults(g_fLog);

	CloseStdioFile(Peaker::m_fTsv);
	}
