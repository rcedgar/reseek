#pragma once

#include "chaq.h"

static const float BAD_SCORE = -9999;
static const float MIN_SANE_SCORE = -1000;
static const float MAX_SANE_SCORE = 1000;

// active alphabets
// alphabet = feature + alphabet_size + logodds...
class flat_alphas
	{
private:
	flat_alphas();

public:
	static uint32 m_nfeat;
	static uint32 m_entropyfi;
	static vector<string> m_feature_names;
	static uint32_t *m_alpha_sizes;
	static float **m_unweighted_logoddsvec;
	static float **m_weighted_logoddsvec;
	static float *m_weights;
	static uint32_t *m_feature_block_offsets;
	static uint32_t m_sum_alpha_sizes;
	static uint32_t m_compound_alpha_size;
	static uint32_t *m_axes;
	static vector<string> m_symbolsvec;

public:
	static uint get_nfeat() { assert(m_nfeat != 0); return m_nfeat; }
	
	static void alloc(uint32 nfeat);

	static void init(const vector<string> &feature_names);

	static void read_logoddsvec(const vector<string> &fns);

	static void read_logoddsvec_pattern(const string &fnpattern,
		const vector<string> &feature_names,
		const vector<uint> &alpha_sizes);

	static void read_logoddsvec_pattern(const string &fnpattern);

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

	static void load_alphas_obsolete(
		const vector<string> &feature_names,
		const string &logoddsfnpattern);

	static void init_from_alphadir(const vector<string> &alpha_names);
	};