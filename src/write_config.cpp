#include "myutils.h"
#include "flat_bench.h"
#include "flat_params.h"
#include "tabbedlines.h"

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values);

void flat_params::read_config(const string &fn)
	{
	asserta(fn != "");
	FILE *f = OpenStdioFile(fn);
	unset_all_params();

	vector<string> lines;
	ReadLinesFromFile(fn, lines);
	tabbedlines tl(lines);

#define x(name)	m_##name = tl.get_float(#name);
	x(open);
	x(ext);
	x(self_w);
	x(rev_w);
	x(lddt_w);
	x(lddtx_w);
	x(lddtpow_w);
	x(dali_w);
	x(entropy_w);
	x(rotfreetm_w);
	x(LDDT_R0);
#undef x

#define x(name)	m_##name = tl.get_int(#name);
	x(nn_min_offset);
	x(distmx_bandwidth);
	x(turnd_w);
	x(angle_n);
	x(LDDT_nr_thresholds);
#undef x

	m_LDDT_thresholds =
		tl.get_float_vec("LDDT_thresholds", m_LDDT_nr_thresholds);

	m_nfeat = tl.get_int("nfeat");
	tl.get_str_vec("alpha_names", flat_params::m_feature_names);

	uint n;
	flat_params::m_alpha_sizes = tl.get_int_vec("alpha_sizes", n);
	asserta(n == m_nfeat);

	flat_params::m_undef_values = tl.get_int16_vec("undef_values", n);
	asserta(n == m_nfeat);

	flat_params::m_undef_codes = tl.get_int8_vec("undef_codes", n);
	asserta(n == m_nfeat);

	flat_params::m_weights = tl.get_float_vec("alpha_weights", n);
	asserta(n == m_nfeat);

	asserta(flat_params::m_feature_names.size() == m_nfeat);
	m_unweighted_logoddsvec = myalloc(float *, m_nfeat);
	m_weighted_logoddsvec = myalloc(float *, m_nfeat);
	m_thresholds = myalloc(uint16_t *, m_nfeat);
	for (uint fi = 0; fi < m_nfeat; ++fi)
		{
		const string &feature_name = flat_params::m_feature_names[fi];
		string fan_name;
		get_fan_name(feature_name, fan_name);
		FAN fan = str2FAN(fan_name.c_str());
		m_fans.push_back(fan);

		uint AS = m_alpha_sizes[fi];
		const string name = "unweighted_logodds_" + feature_name;
		uint alpha_size = m_alpha_sizes[fi];
		uint AS2;
		m_unweighted_logoddsvec[fi] = tl.get_float_flat_square_mx(name, AS2);
		m_weighted_logoddsvec[fi] = myalloc(float, AS*AS);
		asserta(AS2 == alpha_size);
		if (flat_params::feature_is_binned(fan))
			{
			const string tname = "thresholds_" + feature_name;
			m_thresholds[fi] = tl.get_int16_flat_vec(tname, alpha_size-1);
			}
		else
			m_thresholds[fi] = 0;
		}

	CloseStdioFile(f);
	}

void flat_params::post_config_setup()
	{
	m_sum_alpha_sizes = 0;
	m_compound_alpha_size = 1;
	m_entropyfi = UINT_MAX;
	m_axes = myalloc(uint32_t, m_nfeat);//TODO consolidate alloc
	m_feature_block_offsets = myalloc(uint32_t, m_nfeat);

	for (uint fi = 0; fi < m_nfeat; ++fi)
		{
		const string &feature_name = m_feature_names[fi];
		string fan_name;
		get_fan_name(feature_name, fan_name);
		FAN fan = str2FAN(fan_name.c_str());
 		if (StartsWith(feature_name, "sec") || feature_name == "Conf")
			m_entropyfi = fi;
		uint alpha_size = m_alpha_sizes[fi];
		m_sum_alpha_sizes += alpha_size;
		m_axes[fi] = m_compound_alpha_size;
		m_compound_alpha_size *= alpha_size;
		}

	set_feature_block_offsets();
	}

void flat_params::load_config(const string &fn)
	{
	read_config(fn);
	post_config_setup();
	apply_current_weights();
	}

void flat_params::write_config(const string &fn)
	{
	if (fn == "") return;
	const uint nfeat = flat_params::get_nfeat();
	FILE *f = CreateStdioFile(fn);

	tabbedlines tl;
#define x(name)	tl.put_float(#name, m_##name)
	x(open);
	x(ext);
	x(self_w);
	x(rev_w);
	x(lddt_w);
	x(lddtx_w);
	x(lddtpow_w);
	x(dali_w);
	x(entropy_w);
	x(rotfreetm_w);
	x(LDDT_R0);
#undef x

#define x(name)	tl.put_int(#name, m_##name)
	x(nn_min_offset);
	x(distmx_bandwidth);
	x(turnd_w);
	x(angle_n);
	x(LDDT_nr_thresholds);
#undef x

	tl.put_float_vec("LDDT_thresholds",
		m_LDDT_thresholds, m_LDDT_nr_thresholds);

	tl.put_int("nfeat", nfeat);
	tl.put_str_vec("alpha_names", flat_params::m_feature_names);
	tl.put_int_vec("alpha_sizes", flat_params::m_alpha_sizes, nfeat);
	tl.put_int16_vec("undef_values", flat_params::m_undef_values, nfeat);
	tl.put_int8_vec("undef_codes", flat_params::m_undef_codes, nfeat);
	tl.put_float_vec("alpha_weights", flat_params::m_weights, nfeat);

	for (uint fi = 0; fi < nfeat; ++fi)
		{
		const string &feature_name = flat_params::m_feature_names[fi];
		string fan_name;
		get_fan_name(feature_name, fan_name);
		FAN fan = str2FAN(fan_name.c_str());

		uint AS = m_alpha_sizes[fi];
		const string name = "unweighted_logodds_" + feature_name;
		tl.put_float_flat_square_mx(name, AS, m_unweighted_logoddsvec[fi]);
		if (flat_params::feature_is_binned(fan))
			{
			const string tname = "thresholds_" + feature_name;
			cp_uint16_t thresholds =
				chaq::get_hard_coded_thresholds(fan, AS);
			tl.put_int16_flat_vec(tname, thresholds, AS-1);
			}
		}

	tl.to_tsv(f);
	CloseStdioFile(f);
	}

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

	chaq::m_enable_hard_coded_parameters = true;

	vector<string> feature_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		feature_names, weights,
		scalar_names, scalar_values);//TODO should be in flat_params

	const size_t n = feature_names.size();
	asserta(weights.size() == n);
	unordered_map<string, float> name2weight;
	for (size_t i = 0; i < n; ++i)
		name2weight[feature_names[i]] = weights[i];

	flat_params::set_params(scalar_names, scalar_values);
	flat_params::load_alphas(feature_names, opt(mxpattern));
	flat_params::apply_weights(name2weight);

	flat_params::write_config(opt(output));
	flat_params::read_config(opt(output));
	flat_params::write_config(opt(output2));
	}
