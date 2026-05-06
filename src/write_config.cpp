#include "myutils.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "flat_bench.h"

void cmd_write_config()
	{
	asserta(optset_mxpattern);
	asserta(optset_output);

	asserta(!optset_lookup);
	asserta(!optset_fapattern);
	asserta(!optset_spec);
	asserta(!optset_varstr);

	const string &VarStr = g_Arg1;

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

	const size_t n = feature_names.size();
	asserta(weights.size() == n);
	unordered_map<string, float> name2weight;
	for (size_t i = 0; i < n; ++i)
		name2weight[feature_names[i]] = weights[i];

	flat_features::load_alphas(feature_names, opt(mxpattern));
	flat_features::apply_weights(name2weight);

	FILE *f = CreateStdioFile(opt(output));
	flat_features::write_config(f);
	flat_params::write_config(f);
	CloseStdioFile(f);
	}
