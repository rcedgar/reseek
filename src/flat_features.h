#pragma once

#include "chaq.h"

static const float BAD_SCORE = -9999;
static const float MIN_SANE_SCORE = -1000;
static const float MAX_SANE_SCORE = 1000;

class flat_features
	{
public:
	uint32 m_nfeat = 0;
	vector<string> m_feature_names;
	uint32_t *m_alpha_sizes = 0;
	float **m_unweighted_logoddsvec = 0;
	float **m_weighted_logoddsvec = 0;
	float *m_weights = 0;
	uint32_t *m_feature_block_offsets = 0;
	uint32_t m_sum_alpha_sizes = 0;
	uint32_t m_compound_alpha_size = 0;
	uint32_t *m_axes = 0;
	vector<string> m_symbolsvec;

public:
	uint get_nfeat() const { asserta(m_nfeat != 0); return m_nfeat; }
	
	void alloc(uint32 nfeat);

	void init(const vector<string> &feature_names);

	void read_logoddsvec(const vector<string> &fns);

	void read_logoddsvec_pattern(const string &fnpattern,
		const vector<string> &feature_names,
		const vector<uint> &alpha_sizes);

	void read_logoddsvec_pattern(const string &fnpattern);

	void apply_weights(const vector<float> &weights);

	void apply_weights(const unordered_map<string, float> &name2weight);

	void apply_unit_weights();

	const string &get_symbols(uint fi);

	float prof_col_score(
		const uint8_t *profQ, uint LQ, uint posQ,
		const uint8_t *profT, uint LT, uint posT) const;

	const uint32_t *get_feature_block_offsets() const;

	uint32_t get_compound_alpha_size() const { return m_compound_alpha_size; }

	uint32_t get_sum_alpha_sizes() const
		{
		assert(m_sum_alpha_sizes > 0);
		return m_sum_alpha_sizes;
		}
	
	const unsigned char *get_letter2char(uint fi) const
		{
		assert(fi < m_nfeat);
		uint alpha_size = m_alpha_sizes[fi];
		return chaq::get_letter2char(alpha_size);
		}

	const uint8_t *get_char2letter(uint fi) const
		{
		assert(fi < m_nfeat);
		uint alpha_size = m_alpha_sizes[fi];
		return chaq::get_char2letter(alpha_size);
		}

	uint8_t component_codes_to_compound_code(
		const vector<uint8_t> &codes) const;

	void compound_code_to_component_codes(
		uint8_t code, vector<uint8_t> &codes) const;

	float get_compound_subst_score_slow(
		uint8_t code1, uint8_t code2) const;

	void check_sane_scores() const;

	void set_symbolsvec();
	void set_feature_block_offsets();

public:
	static void get_logodds_symbols(const float *logodds,
		uint alpha_size, string &symbols);

	static uint read_logodds(const string &fn, vector<float> &logodds);
	static uint lines2logoddsmx(const vector<string> &lines,
		vector<float> &logoddsmx);

	static void make_logoddsfn_pattern(
		const string &fnpattern,
		const string &feature_name,
		uint alpha_size,
		string &fn);
	};