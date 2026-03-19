#include "myutils.h"
#include "tabbedlines.h"
#include "flat_helpers.h"
#include "flat_features.h"

void flat_features::init(
	const vector<string> &feature_names,
	const vector<uint> &alpha_sizes)
	{
	asserta(m_nfeat == 0);
	alloc(uint(feature_names.size()));
	asserta(alpha_sizes.size() == m_nfeat);
	m_feature_names = feature_names;
	memcpy(m_alpha_sizes, alpha_sizes.data(),
		m_nfeat*sizeof(m_alpha_sizes[0]));
	}

uint flat_features::lines2logoddsmx(
	const vector<string> &lines,
	vector<float> &logoddsmx)
	{
	logoddsmx.clear();
	tabbedlines tl(lines);
	uint alpha_size = tl.get_int("logodds");
	asserta(alpha_size != 0);
	logoddsmx.resize(alpha_size*alpha_size);
	tl.get_float_flat_square_mx(alpha_size, logoddsmx.data());
	return alpha_size;
	}

uint flat_features::read_logodds(
	const string &fn,
	vector<float> &logoddsmx)
	{
	vector<string> lines;
	ReadLinesFromFile(fn, lines);
	return lines2logoddsmx(lines, logoddsmx);
	}

void flat_features::alloc(uint32 nfeat)
	{
	if (nfeat == m_nfeat)
		return;
	assert(nfeat > 0);
	assert(m_nfeat == 0);
	assert(m_weights == 0);
	assert(m_unweighted_logoddsvec == 0);
	assert(m_weighted_logoddsvec == 0);
	assert(m_feature_block_offsets == 0);

	m_nfeat = nfeat;
	m_weights = myalloc(float, m_nfeat);
	m_alpha_sizes = myalloc(uint32_t, m_nfeat);
	m_unweighted_logoddsvec = myalloc(float *, m_nfeat);
	m_weighted_logoddsvec = myalloc(float *, m_nfeat);
	m_feature_block_offsets = myalloc(uint32_t, m_nfeat);
	}

void flat_features::read_logoddsvec(const vector<string> &fns)
	{
	_chkmem();//@@
	uint nfeat = uint(fns.size());
	alloc(nfeat);
	_chkmem();//@@
	for (uint i = 0; i < m_nfeat; ++i)
		{
		const string &fn = fns[i];
		vector<float> logodds;
		uint AS = read_logodds(fn, logodds);
		m_alpha_sizes[i] = AS;
		uint bytes = AS*AS*sizeof(float);
		m_unweighted_logoddsvec[i] = myalloc(float, bytes);
		m_weighted_logoddsvec[i] = myalloc(float, bytes);
		for (uint k = 0; k < AS*AS; ++k)
			{
			float score = logodds[k];
			assert(score >= MIN_SANE_SCORE && score <= MAX_SANE_SCORE);
			m_unweighted_logoddsvec[i][k] = score;
			m_weighted_logoddsvec[i][k] = BAD_SCORE;
			}
		_chkmem();//@@
		}
	}

void flat_features::read_logoddsvec_pattern(
	const string &fnpattern,
	const vector<string> &feature_names,
	const vector<uint> &alpha_sizes)
	{
	uint nfeat = SIZE(feature_names);
	asserta(SIZE(alpha_sizes) == nfeat);
	alloc(nfeat);

	m_feature_names = feature_names;

	memcpy(m_alpha_sizes, alpha_sizes.data(),
		nfeat*sizeof(m_alpha_sizes[0]));

	vector<string> fns(m_nfeat);
	for (uint fi = 0; fi < m_nfeat; ++fi)
		make_logoddsfn_pattern(
			fnpattern,
			feature_names[fi],
			alpha_sizes[fi],
			fns[fi]);
	read_logoddsvec(fns);
	}

void flat_features::read_logoddsvec_pattern(
	const string &fnpattern)
	{
	asserta(m_nfeat > 0);
	vector<string> fns(m_nfeat);
	for (uint fi = 0; fi < m_nfeat; ++fi)
		{
		_chkmem();//@@
		make_logoddsfn_pattern(
			fnpattern,
			m_feature_names[fi],
			m_alpha_sizes[fi],
			fns[fi]);
		_chkmem();//@@
		}
	read_logoddsvec(fns);
	_chkmem();//@@
	}

void flat_features::check_sane_scores() const
	{
	for (uint fi = 0; fi < m_nfeat; ++fi)
		{
		uint AS = m_alpha_sizes[fi];
		const float *low = m_weighted_logoddsvec[fi];
		const float *lou = m_unweighted_logoddsvec[fi];
		for (uint k = 0; k < AS*AS; ++k)
			{
			float wscore = low[k];
			float uscore = lou[k];
			asserta(wscore >= MIN_SANE_SCORE && wscore <= MAX_SANE_SCORE);
			asserta(uscore >= MIN_SANE_SCORE && uscore <= MAX_SANE_SCORE);
			}
		}
	}

void flat_features::set_symbolsvec()
	{
	asserta(m_nfeat > 0);
	m_symbolsvec.clear();
	m_symbolsvec.resize(m_nfeat);
	for (uint fi = 0; fi < m_nfeat; ++fi)
		get_logodds_symbols(
			m_unweighted_logoddsvec[fi],
			m_alpha_sizes[fi],
			m_symbolsvec[fi]);
	}

const string &flat_features::get_symbols(uint fi)
	{
	asserta(fi < m_nfeat);
	if (m_symbolsvec.empty())
		set_symbolsvec();
	asserta(m_symbolsvec.size() == m_nfeat);
	return m_symbolsvec[fi];
	}

// @=name, %=AS
void flat_features::make_logoddsfn_pattern(
	const string &fnpattern,
	const string &feature_name,
	uint alpha_size,
	string &fn)
	{
	fn.clear();
	for (auto c : fnpattern)
		{
		if (c == '@')
			fn += feature_name;
		else if (c == '%')
			fn += to_string(alpha_size);
		else
			fn += c;
		}
	}

void flat_features::get_logodds_symbols(
	const float *logodds, uint alpha_size, string &symbols)
	{
	symbols.clear();
	float min_score = FLT_MAX;
	float max_score = FLT_MAX;
	for (uint i = 0; i < alpha_size*alpha_size; ++i)
		{
		float score = logodds[i];
		min_score = (i == 0 ? score : min(min_score, score));
		max_score = (i == 0 ? score : max(max_score, score));
		}

	// __. +*^
	// 0123456
	static const char s[7] = { 'V', '_', '.', ' ', '+', '*', '^' };
	for (uint i = 0; i < alpha_size; ++i)
		{
		for (uint j = 0; j < alpha_size; ++j)
			{
			float score = logodds[alpha_size*i + j];
			assert(score >= min_score && score <= max_score);
			uint k = uint(7*(score - min_score)/(max_score - min_score + max_score/7));
			symbols += s[k];
			}
		}
	}

void flat_features::set_feature_block_offsets()
	{
	assert(m_feature_block_offsets != 0);
	assert(m_nfeat > 0);
	m_sum_alpha_sizes = 
		get_flat_pssm_feature_block_offsets(
			m_nfeat, m_alpha_sizes, m_feature_block_offsets);
	}

const uint32_t *flat_features::get_feature_block_offsets() const
	{
	assert(m_feature_block_offsets != 0);
	return m_feature_block_offsets;
	}

void flat_features::apply_weights(const vector<float> &weights)
	{
	assert(m_weights != 0);
	memcpy(m_weights, weights.data(), m_nfeat*sizeof(float));
	float sumw = 0;
	for (uint i = 0; i < m_nfeat; ++i) sumw += m_weights[i];
	asserta(sumw > 1e-6);
	for (uint i = 0; i < m_nfeat; ++i) m_weights[i] /= sumw;

	for (uint fi = 0; fi < m_nfeat; ++fi)
		{
		uint AS = m_alpha_sizes[fi];
		uint N = AS*AS;
		for (uint k = 0; k < N; ++k)
			{
			float uwscore = m_unweighted_logoddsvec[fi][k];
			assert(uwscore >= MIN_SANE_SCORE && uwscore <= MAX_SANE_SCORE);

			float wscore = uwscore*m_weights[fi];
			assert(wscore >= MIN_SANE_SCORE && wscore <= MAX_SANE_SCORE);

			m_weighted_logoddsvec[fi][k] = wscore;
			}
		}
	check_sane_scores();
	}

void flat_features::apply_unit_weights()
	{
	vector<float> w(m_nfeat, 1);
	apply_weights(w);
	}

float flat_features::prof_col_score(
	const uint8_t *profQ, uint LQ, uint posQ,
	const uint8_t *profT, uint LT, uint posT) const
	{
	assert(posQ < LQ);
	assert(posT < LT);
	float score = 0;
	for (uint32_t fi = 0; fi < m_nfeat; ++fi)
		{
		const uint32_t AS_fi = m_alpha_sizes[fi];
		const float *logodds_fi = m_weighted_logoddsvec[fi];
		const uint8_t *profQ_fi = profQ + fi*LQ;
		const uint8_t *profT_fi = profT + fi*LT;
		const uint8_t codeQ = profQ_fi[posQ];
		const uint8_t codeT = profT_fi[posT];
		const float *logodds_row = logodds_fi + codeQ*AS_fi;
		score += logodds_row[codeT];
		}
	return score;
	}
