#pragma once

#include "chaq.h"

static const float BAD_SCORE = -9999;
static const float MIN_SANE_SCORE = -1000;
static const float MAX_SANE_SCORE = 1000;

class flat_params
	{
private:
	flat_params();

public:
	static uint32 m_nfeat;
	static uint32 m_entropyfi;
	static vector<string> m_feature_names;
	static vector<FAN> m_fans;
	static uint32_t *m_alpha_sizes;
	static uint16_t *m_undef_values;
	static uint8_t *m_undef_codes;
	static uint16_t **m_thresholds;
	static float **m_unweighted_logoddsvec;
	static float **m_weighted_logoddsvec;
	static float *m_weights;
	static uint32_t *m_feature_block_offsets;
	static uint32_t m_sum_alpha_sizes;
	static uint32_t m_compound_alpha_size;
	static uint32_t *m_axes;
	static vector<string> m_symbolsvec;

	static float m_open;
	static float m_ext;
	static float m_self_w;
	static float m_rev_w;
	static float m_lddt_w;
	static float m_lddtx_w;
	static float m_lddtpow_w;
	static float m_dali_w;
	static float m_dalix_w;
	static float m_entropy_w;
	static float m_rotfreetm_w;

	static uint32_t m_nn_min_offset;
	static uint32_t m_distmx_bandwidth;
	static uint32_t m_turnd_w;
	static uint32_t m_angle_n;

	static float m_LDDT_R0;
	static const float *m_LDDT_thresholds;
	static uint m_LDDT_nr_thresholds;

public:
	static void unset_all_params()
		{
		for (uint fi = 0; fi < m_nfeat; ++fi)
			{
			myfree(m_unweighted_logoddsvec[fi]);
			myfree(m_weighted_logoddsvec[fi]);
			myfree(m_thresholds[fi]);
			}

		myfree(m_unweighted_logoddsvec);
		myfree(m_weighted_logoddsvec);
		myfree(m_thresholds);

		m_unweighted_logoddsvec = 0;
		m_weighted_logoddsvec = 0;
		m_thresholds = 0;

		m_open = 0;
		m_ext = 0;
		m_self_w = 0;
		m_rev_w = 0;
		m_lddt_w = 0;
		m_lddtx_w = 0;
		m_lddtpow_w = 0;
		m_dali_w = 0;
		m_dalix_w = 0;
		m_entropy_w = 0;
		m_rotfreetm_w = 0;

		m_nfeat = 0;
		m_entropyfi = UINT_MAX;
		m_sum_alpha_sizes = 0;
		m_sum_alpha_sizes = 0;
		m_compound_alpha_size = 0;
		m_nn_min_offset = 0;
		m_distmx_bandwidth = 0;
		m_turnd_w = 0;
		m_angle_n = 0;

		m_feature_names.clear();
		m_fans.clear();
		m_symbolsvec.clear();

		m_LDDT_R0 = 0;
		m_LDDT_thresholds = 0;
		m_LDDT_nr_thresholds = 0;

#define x(name)	if (m_##name != 0) { myfree(m_##name); m_##name = 0; }
		x(alpha_sizes);
		x(undef_values);
		x(feature_block_offsets);
		x(axes);
		x(weights);
#undef x
		}

	static uint get_nfeat() { assert(m_nfeat != 0); return m_nfeat; }
	
	static void alloc(uint32 nfeat);

	static void init(const vector<string> &feature_names);

	static void read_logoddsvec(const vector<string> &fns);

	static void read_logoddsvec_pattern(const string &fnpattern,
		const vector<string> &feature_names,
		const vector<uint> &alpha_sizes);

	static void read_logoddsvec_pattern(const string &fnpattern);

	static void apply_current_weights();
	static void apply_weights(const vector<float> &weights);

	static void apply_weights(const unordered_map<string, float> &name2weight);

	static void apply_unit_weights();

	static const string &get_symbols(uint fi);

	static float prof_col_score(
		const uint8_t *profQ, uint LQ, uint posQ,
		const uint8_t *profT, uint LT, uint posT);

	static const uint32_t *get_feature_block_offsets();

	static uint32_t get_compound_alpha_size() { return m_compound_alpha_size; }

	static uint32_t get_sum_alpha_sizes()
		{
		assert(m_sum_alpha_sizes > 0);
		return m_sum_alpha_sizes;
		}
	
	static const unsigned char *get_letter2char(uint fi)
		{
		assert(fi < m_nfeat);
		uint alpha_size = m_alpha_sizes[fi];
		return chaq::get_letter2char(alpha_size);
		}

	static const uint8_t *get_char2letter(uint fi)
		{
		assert(fi < m_nfeat);
		uint alpha_size = m_alpha_sizes[fi];
		return chaq::get_char2letter(alpha_size);
		}

	static uint8_t component_codes_to_compound_code(
		const vector<uint8_t> &codes);

	static void compound_code_to_component_codes(
		uint8_t code, vector<uint8_t> &codes);

	static float get_compound_subst_score_slow(
		uint8_t code1, uint8_t code2);

	static void check_sane_scores();

	static void set_symbolsvec();
	static void set_feature_block_offsets();

	static void get_compound_logodds_slow(vector<float> &logodds);

	static void get_logodds_symbols(const float *logodds,
		uint alpha_size, string &symbols);

	static uint read_logodds(const string &fn, vector<float> &logodds);
	static void write_logodds(const string &fn,
		const vector<float> &logodds, uint alpha_size);
	static void logodds2lines(const vector<float> &logodds,
		uint alpha_size, vector<string> &lines);
	static uint lines2logoddsmx(const vector<string> &lines,
		vector<float> &logoddsmx);

	static void make_logoddsfn_pattern(
		const string &fnpattern,
		const string &feature_name,
		uint alpha_size,
		string &fn);

	static void load_alphas(
		const vector<string> &feature_names,
		const string &logoddsfnpattern);

	//static void set_alphas(const vector<string> &feature_names);
	static void write_config(const string &fn);
	static void read_config(const string &fn);
	static void load_config(const string &fn);
	static void post_config_setup();

	static void get_fan_name(
		const string &feature_name, 
		string &fan_name);
	static void set_params(
		const vector<string> &names,
		const vector<float> &values);

	static bool need_reverse();
	static bool need_distmx();
	static bool need_prof();
	static bool need_self();
	static bool need_alignx();

	static FAN flat_params::get_fan(uint fi)
		{
		assert(fi < m_nfeat);
		return m_fans[fi];
		}

	static uint flat_params::get_alpha_size(uint fi)
		{
		assert(fi < m_nfeat);
		return m_alpha_sizes[fi];
		}

	static uint16_t flat_params::get_undef_value(uint fi)
		{
		assert(fi < m_nfeat);
		return m_undef_values[fi];
		}

	static uint8_t flat_params::get_undef_code(uint fi)
		{
		assert(fi < m_nfeat);
		return m_undef_codes[fi];
		}

	static cp_uint16_t get_thresholds(uint fi)
		{
		assert(fi < m_nfeat);
		return m_thresholds[fi];
		}

	static bool feature_is_binned_fi(uint fi)
		{
		assert(fi < m_nfeat);
		FAN fan = m_fans[fi];
		return feature_is_binned(fan);
		}

	static bool feature_is_binned(FAN fan)
		{
		switch (fan)
			{
		case FAN_nendist:
		case FAN_rendist:
		case FAN_pendist:
		case FAN_mendist:
		case FAN_fendist:
		case FAN_pack:
		case FAN_ppack:
		case FAN_mpack:
		case FAN_angle:
		case FAN_pmdd:
		case FAN_pmdiff:
		case FAN_turnd:
			return true;
			}
		return false;
		}
	};