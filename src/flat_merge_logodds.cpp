#include "myutils.h"
#include "seqdb.h"
#include "flat_helpers.h"
#include "flat_bench.h"

// Make logodds for compound alphabet using
//   varstr to support non-uniform weights
void cmd_flat_merge_logodds()
	{
	asserta(optset_mxpattern);

	asserta(!optset_lookup);
	asserta(!optset_fapattern);
	asserta(!optset_varstr);

	const string &VarStr = g_Arg1;

	flat_bench FB;
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

	flat_features ff;
	ff.init(feature_names);
	ff.read_logoddsvec_pattern(opt(mxpattern));
	ff.apply_weights(weights);
	const uint compound_alpha_size = ff.get_compound_alpha_size();
	ProgressLog("compound alpha_size %u\n", compound_alpha_size);

	vector<float> logodds;
	ff.get_compound_logodds_slow(logodds);
	FILE *f = CreateStdioFile(opt(output));
	write_flat_logoddsmx(f, logodds, compound_alpha_size, false);
	CloseStdioFile(f);
	ProgressLog("Ok\n");
	}
