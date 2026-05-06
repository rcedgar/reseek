#include "myutils.h"
#include "flat_params.h"
#include "flat_bench.h"
#include "flat_helpers.h"
#include "flat_bench_struct_feature.h"

float flat_rotfreetm(
	uint loQ, uint LQ,
	uint loT, uint LT,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	const string &path)
	{
	uint nmatch = 0;
	for (auto c : path) if (c == 'M') ++nmatch;
	const uint M = flat_params::m_distmx_bandwidth;

	vector<uint32_t> posQs;
	vector<uint32_t> posTs;
	path2posvecs(path, loQ, LQ, loT, LT, posQs, posTs);
	const uint K = uint(posQs.size());
	Die("TODO");
	return 0;
	}
	
float flat_bench_struct_feature::get_rotfreetm(uint idxQ, uint idxT,
	const flat_aligner &fa) const
	{
	asserta(idxQ < m_distmxs.size());
	asserta(idxT < m_distmxs.size());

	const sid_t *distmxQ = m_distmxs[idxQ];
	const sid_t *distmxT = m_distmxs[idxT];

	string path;
	fa.get_path_str(path);
	float score = flat_rotfreetm(
		fa.m_loQ, fa.m_LQ, fa.m_loT, fa.m_LT,
		distmxQ, distmxT, path);
	return score;
	}

void flat_bench::align_pair_rotfreetm(flat_aligner &fa,
	const string &labelQ, const string &labelT)
	{
	fa.alloc();

	string labQ, labT;
	trunc_label(labelQ, labQ);
	trunc_label(labelT, labT);
	uint DomIdxQ = UINT_MAX;
	uint DomIdxT = UINT_MAX;
	for (uint i = 0; i < uint(m_Labels.size()); ++i)
		{
		if (m_Labels[i] == labQ)
			DomIdxQ = i;
		if (m_Labels[i] == labT)
			DomIdxT = i;
		}
	asserta(DomIdxQ != UINT_MAX);
	asserta(DomIdxT != UINT_MAX);

	const uint8_t *profT = m_fp.get_profile(DomIdxT);
	const uint LT = m_fp.get_length(DomIdxT);
	fa.cacheT(labelT, profT, LT);

	const uint8_t *profQ = m_fp.get_profile(DomIdxQ);
	const uint LQ = m_fp.get_length(DomIdxQ);
	fa.alignQ(labelQ, profQ, LQ);
	float score = fa.m_score;
	asserta(!flat_params::need_alignx());
	fa.write_aln(g_fLog);

	const sid_t *distmxQ = m_distmxs[DomIdxQ];
	const sid_t *distmxT = m_distmxs[DomIdxT];

	string path;
	fa.get_path_str(path);
	float rotfreetm = flat_rotfreetm(
		fa.m_loQ, fa.m_LQ, fa.m_loT, fa.m_LT,
		distmxQ, distmxT, path);
	Log("rotfreetm=%.4f\n", rotfreetm);
	}

void cmd_rotfreetm()
	{
	asserta(optset_lookup);
	asserta(optset_fapattern);
	asserta(optset_mxpattern);
	asserta(optset_input);
	asserta(optset_label1);
	asserta(optset_label2);

	asserta(!optset_spec);
	asserta(!optset_varstr);

	const string &VarStr = g_Arg1;

	flat_bench FB;
	FB.ReadLookup(opt(lookup));
	if (optset_dope)
		FB.ReadDope(opt(dope));

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

	flat_params::load_alphas(feature_names, opt(mxpattern));

	vector<flat_chain_t *> chains;
	read_flat_chains(opt(input), chains);

	FB.load_profiles(opt(fapattern));
	FB.set_distmxs(chains);
	FB.UpdateParamsFromVarStr(VarStr);

	flat_aligner fa;
	FB.align_pair_rotfreetm(fa, opt(label1), opt(label2));
	}