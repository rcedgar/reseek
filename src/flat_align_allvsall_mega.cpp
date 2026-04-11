#include "myutils.h"
#include "flat_bench.h"
#include "alpha.h"
#include "sort.h"
#include "triangle.h"
#include "flat_helpers.h"
#include "flat_aligner.h"
#include <unordered_map>
#include <unordered_set>

#define	SHOW_PROGRESS	1

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values);

void cmd_flat_align_allvsall_mega()
	{
	asserta(optset_lookup);
	asserta(optset_fapattern);
	asserta(optset_mxpattern);

	asserta(!optset_dope);
	asserta(!optset_spec);
	asserta(!optset_varstr);

	const string &VarStr = g_Arg1;

	flat_bench FB;
	FB.ReadLookup(opt(lookup));

	vector<string> param_names;
	vector<float> param_values;
	ParseVarStr(VarStr, param_names, param_values);

	vector<string> feature_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		feature_names, weights,
		scalar_names, scalar_values);

	FB.load_alphas_and_profiles(
		feature_names, weights, opt(fapattern), opt(mxpattern));
	FB.UpdateParamsFromVarStr(VarStr);
	FB.ProgressLogParams();
	FB.Alloc();

	FB.m_f_tsv_all_vs_all = CreateStdioFile(opt(output));
	FB.Search("allvsall");
	FB.SetScoreOrder();
	FB.Bench();
	//FB.WriteHits(opt(output), true);
	}
