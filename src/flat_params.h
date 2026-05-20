#pragma once

#include "chaq.h"
#include "collect.h"
#include "fan.h"

static const float BAD_SCORE = -9999;
static const float MIN_SANE_SCORE = -1000;
static const float MAX_SANE_SCORE = 1000;

class flat_params
	{
public:
	// static fixed parameters
	// changing these requires re-training
	// log-odds and quantization thresholds
	static uint32_t m_nn_min_offset;
	static uint32_t m_distmx_bandwidth;
	static uint32_t m_turnd_w;
	static uint32_t m_angle_n;

	static float m_LDDT_R0;
	static const float *m_LDDT_thresholds;
	static uint m_LDDT_nr_thresholds;

public:
	// alignment
	static float m_open;
	static float m_ext;

	// test statistic
	static float m_self_w;
	static float m_rev_w;
	static float m_lddt_w;
	static float m_lddtx_w;
	static float m_dali_w;
	static float m_dalix_w;
	static float m_nurev_w;

	// filters
	static float m_mega_filter_min_fwd;
	static float m_nu_filter_self_w;
	static float m_nu_filter_rev_w;
	static float m_nu_filter_min_fwd_score;
	static float m_nu_filter_min_combined_score;

	// alphabets
	static uint32 m_nfeat;
	static vector<string> m_alpha_names;
	static vector<FAN> m_fans;
	static uint32_t *m_alpha_sizes;
	static float **m_unweighted_logoddsvec;
	static float **m_weighted_logoddsvec;
	static float *m_weights;
	static uint32_t *m_feature_block_offsets;
	static uint32_t m_sum_alpha_sizes;
	static uint32_t m_compound_alpha_size;
	static uint32_t *m_axes;
	static uint16_t *m_medians;
	static uint16_t **m_thresholds;
	static vector<string> m_symbolsvec;

public:
	void set_params(
		const vector<string> &names,
		const vector<float> &values);

	bool need_reverse();
	bool need_distmx();
	bool need_self();
	bool need_nu_self();
	bool need_alignx();

	void logme();

///////////////////
// from flat_alphas
///////////////////
	void set_alpha_names(const vector<string> &alpha_names);

	uint get_nfeat() { assert(m_nfeat != 0); return m_nfeat; }
	
	void alloc(uint32 nfeat);

	void read_logoddsvec(const vector<string> &fns);

	void apply_weights(const vector<float> &weights);

	void apply_weights(const unordered_map<string, float> &name2weight);

	void apply_unit_weights();

	const string &get_symbols(uint fi);

	float prof_col_score(
		const uint8_t *profQ, uint LQ, uint posQ,
		const uint8_t *profT, uint LT, uint posT) const;

	const uint32_t *get_feature_block_offsets() const;

	uint32_t get_compound_alpha_size() const { return m_compound_alpha_size; }

	uint32_t get_fi(FAN fan, uint alpha_size, bool errok = false) const
		{
		for (uint i = 0; i < m_nfeat; ++i)
			if (m_fans[i] == fan && m_alpha_sizes[i] == alpha_size)
				return i;
		if (!errok)
			Die("get_fi(%u=%s, alpha_size=%u)",
				fan, FAN2str(fan), alpha_size);
		return UINT_MAX;
		}

	uint32_t get_sum_alpha_sizes()
		{
		assert(m_sum_alpha_sizes > 0);
		return m_sum_alpha_sizes;
		}
	
	const unsigned char *get_letter2char(uint fi)
		{
		assert(fi < m_nfeat);
		uint alpha_size = m_alpha_sizes[fi];
		return chaq::get_letter2char(alpha_size);
		}

	const uint8_t *get_char2letter(uint fi)
		{
		assert(fi < m_nfeat);
		uint alpha_size = m_alpha_sizes[fi];
		return chaq::get_char2letter(alpha_size);
		}

	uint8_t component_codes_to_compound_code(
		const vector<uint8_t> &codes);

	void compound_code_to_component_codes(
		uint8_t code, vector<uint8_t> &codes);

	float get_compound_subst_score_slow(
		uint8_t code1, uint8_t code2);

	void check_sane_scores();

	void set_symbolsvec();
	void set_feature_block_offsets();

	void get_compound_logodds_slow(vector<float> &logodds);

	void get_logodds_symbols(const float *logodds,
		uint alpha_size, string &symbols);

	uint read_logodds(const string &fn, vector<float> &logodds);
	void write_logodds(const string &fn,
		const vector<float> &logodds, uint alpha_size);
	void logodds2lines(const vector<float> &logodds,
		uint alpha_size, vector<string> &lines);
	uint lines2logoddsmx(const vector<string> &lines,
		vector<float> &logoddsmx);

	void init_from_alphadir(
		const string &arg_alphadir,
		const vector<string> &alpha_names);

	void init_from_fnprefixes(
		const vector<string> &alpha_names,
		const vector<string> &fnprefixes);

	void init_from_collect(
		const collect &C,
		const vector<string> &alpha_names);
	};