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
	asserta(!optset_scale);	// use -scalef

	float scalef = 1.0f;
	if (optset_scalef)
		scalef = float(opt(scalef));
	const bool as_integers = opt(integers);

	const string &VarStr = g_Arg1;

	flat_bench FB;
	vector<string> param_names;
	vector<float> param_values;
	parse_varstr(VarStr, param_names, param_values);

	vector<string> feature_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	float selfw = 0;
	float revw = 0;
	flat_bench::ClassifyParams(param_names, param_values,
		feature_names, weights,
		scalar_names, scalar_values);

	Die("TODO");
	//flat_params::init(feature_names);
	//flat_params::read_logoddsvec_pattern(opt(mxpattern));
	flat_params::apply_weights(weights);
	const uint compound_alpha_size = flat_params::get_compound_alpha_size();
	const uint AS = compound_alpha_size;
	const uint AS2 = compound_alpha_size*compound_alpha_size;

	ProgressLog("compound alpha_size %u\n", compound_alpha_size);

	vector<float> logodds;
	flat_params::get_compound_logodds_slow(logodds);

	for (uint i = 0; i < AS2; ++i)
		logodds[i] *= scalef;

	if (optset_output)
		{
		ProgressLog("Writing %s\n", opt(output));
		FILE *f = CreateStdioFile(opt(output));
		write_flat_logoddsmx(f, logodds,
			compound_alpha_size, as_integers);
		CloseStdioFile(f);
		}

	if (optset_output2)
		{
		ProgressLog("Writing %s\n", opt(output2));
		FILE *f = CreateStdioFile(opt(output2));
		if (as_integers)
			{
			fprintf(f, "const int logodds[%u] = {\n", AS2);
			for (uint i = 0; i < AS; ++i)
				{
				for (uint j = 0; j < AS; ++j)
					{
					fprintf(f, "%d", int(round(logodds[AS*i + j])));
					if (i+1 < AS || j+1 < AS)
						fprintf(f, ",");
					}
				if (i+1 != AS)
					fprintf(f, "\n");
				}
			fprintf(f, "};\n");
			}
		else
			{
			Die("!integers");
			}
		CloseStdioFile(f);
		}
	}
