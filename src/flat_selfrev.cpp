#include "myutils.h"
#include "flat_chain.h"
#include "flat_bench.h"
#include "flat_params.h"
#include "flat_features.h"
#include "flat_helpers.h"
#include "flat_aligner.h"

void cmd_flat_selfrev()
	{
	asserta(optset_fapattern);
	asserta(optset_mxpattern);
	asserta(optset_output);
	asserta(optset_input);

	asserta(!optset_varstr);

	const string VarStr = g_Arg1;
	const string &chainfn = opt(input);
	vector<string> param_names;
	vector<float> param_values;

	FILE *f = CreateStdioFile(opt(output));

	void ParseVarStr(
		const string &VarStr,
		vector<string> &Names,
		vector<float> &Values);
	ParseVarStr(VarStr, param_names, param_values);

	vector<string> feature_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		feature_names, weights,
		scalar_names, scalar_values);
	flat_params::set_params(scalar_names, scalar_values);

	flat_features::load_alphas(feature_names, opt(mxpattern));
	const uint nfeat = flat_features::get_nfeat();
	unordered_map<string, float> NameToWeight;
	for (uint i = 0; i < nfeat; ++i)
		NameToWeight[feature_names[i]] = weights[i];
	flat_features::apply_weights(NameToWeight);

	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	const uint nchain = uint(chains.size());
	unordered_map<string, uint> label2idx;
	for (uint i = 0; i < nchain; ++i)
		{
		const flat_chain_t *chain = chains[i];
		const string &label = chain->m_label;
		label2idx[label] = i;
		}

	vector<string> fafns(nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		make_fn_pattern(
			opt(fapattern),
			flat_features::m_feature_names[fi],
			fafns[fi]);

	flat_profiles fp;
	fp.read_profiles_from_fastas(fafns, label2idx);

	flat_aligner fa;
	fa.alloc();
	for (uint i = 0; i < nchain; ++i)
		{
		flat_chain_t *chain = chains[i];
		uint L = chain->get_length();
		const uint8_t *prof = fp.get_profile(i);
		float score = fa.get_self_rev_score(chain->m_label, prof, L);
		fprintf(f, "%s\t%.4g\n", chain->m_label.c_str(), score);
		}
	CloseStdioFile(f);
	}