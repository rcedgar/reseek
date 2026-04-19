#include "myutils.h"
#include "flat_bench.h"
#include "alpha.h"
#include "sort.h"
#include "triangle.h"
#include "flat_helpers.h"
#include "flat_aligner.h"
#include "flat_params.h"
#include <unordered_map>
#include <unordered_set>

#define	SHOW_PROGRESS	1

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values);

void cmd_flat_align_selfrev_mega()
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

	flat_features::load_alphas(feature_names, opt(mxpattern));
	FB.load_profiles(opt(fapattern));
	FB.UpdateParamsFromVarStr(VarStr);
	FB.ProgressLogParams();

	FILE *f = CreateStdioFile(opt(output));
	const uint ndom = FB.m_look->get_ndom();
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		ProgressStep(domidx, ndom, "Self-aligning");
		FB.align_pair_selfrev(f, domidx);
		}
	CloseStdioFile(f);
	}
