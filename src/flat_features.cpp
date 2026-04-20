#include "myutils.h"
#include "tabbedlines.h"
#include "flat_helpers.h"
#include "flat_features.h"

uint32 flat_features::m_nfeat;
uint32 flat_features::m_entropyfi = UINT_MAX;
vector<string> flat_features::m_feature_names;
uint32_t *flat_features::m_alpha_sizes;
float **flat_features::m_unweighted_logoddsvec;
float **flat_features::m_weighted_logoddsvec;
float *flat_features::m_weights;
uint32_t *flat_features::m_feature_block_offsets;
uint32_t flat_features::m_sum_alpha_sizes;
uint32_t flat_features::m_compound_alpha_size;
uint32_t *flat_features::m_axes;
vector<string> flat_features::m_symbolsvec;

void flat_features::init(const vector<string> &feature_names)
	{
	asserta(m_nfeat == 0);
	alloc(uint(feature_names.size()));
	m_feature_names = feature_names;
	m_sum_alpha_sizes = 0;
	m_compound_alpha_size = 1;
	m_entropyfi = UINT_MAX;
	for (uint fi = 0; fi < m_nfeat; ++fi)
		{
		const string &feature_name = feature_names[fi];
		if (StartsWith(feature_name, "sec") || feature_name == "Conf")
			m_entropyfi = fi;
		uint alpha_size =
			get_alpha_size_from_feature_name(feature_name);
		m_alpha_sizes[fi] = alpha_size;
		m_sum_alpha_sizes += alpha_size;
		m_axes[fi] = m_compound_alpha_size;
		m_compound_alpha_size *= alpha_size;
		}
	}

void flat_features::logodds2lines(const vector<float> &logodds,
	uint alpha_size, vector<string> &lines)
	{
	tabbedlines tl(lines);
	tl.put_float_flat_square_mx(alpha_size, logodds.data());
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

void flat_features::write_logodds(const string &fn,
		const vector<float> &logodds, uint alpha_size)
	{
	vector<string> lines;
	tabbedlines tl(lines);
	tl.put_float_flat_square_mx(alpha_size, logodds.data());
	tl.to_tsv(fn);
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
	m_axes = myalloc(uint32_t, m_nfeat);
	}

void flat_features::read_logoddsvec(const vector<string> &fns)
	{
	uint nfeat = uint(fns.size());
	alloc(nfeat);
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
		make_logoddsfn_pattern(
			fnpattern,
			m_feature_names[fi],
			m_alpha_sizes[fi],
			fns[fi]);
		}
	read_logoddsvec(fns);
	}

void flat_features::check_sane_scores()
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

const uint32_t *flat_features::get_feature_block_offsets()
	{
	assert(m_feature_block_offsets != 0);
	return m_feature_block_offsets;
	}

void flat_features::apply_weights(
	const unordered_map<string, float> &NameToWeight)
	{
	asserta(SIZE(NameToWeight) == m_nfeat);
	unordered_map<string, uint> NameToIdx;
	for (uint idx = 0; idx < m_nfeat; ++idx)
		NameToIdx[m_feature_names[idx]] = idx;

	for (unordered_map<string, float>::const_iterator iter = NameToWeight.begin();
		iter != NameToWeight.end(); ++iter)
		{
		const string &Name = iter->first;
		float Weight = iter->second;
		unordered_map<string, uint>::const_iterator iter2 =
			NameToIdx.find(Name);
		asserta(iter2 != NameToIdx.end());
		uint idx = iter2->second;
		m_weights[idx] = Weight;

		uint AS = m_alpha_sizes[idx];
		for (uint code = 0; code < AS*AS; ++code)
			m_weighted_logoddsvec[idx][code] =
				m_unweighted_logoddsvec[idx][code]*Weight;
		}
	check_sane_scores();
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
	const uint8_t *profT, uint LT, uint posT)
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

uint8_t flat_features::component_codes_to_compound_code(
	const vector<uint8_t> &component_codes)
	{
	uint compound_code = 0;
	asserta(SIZE(component_codes) == m_nfeat);
	for (uint fi = 0; fi < m_nfeat; ++fi)
		{
		byte component_code = component_codes[fi];
		compound_code += component_code*m_axes[fi];
		}
	byte b = byte(compound_code);
	asserta(uint(b) == compound_code);
	return b;
	}

void flat_features::compound_code_to_component_codes(
	uint8_t compound_code, vector<uint8_t> &codes)
	{
	codes.clear();
	codes.resize(m_nfeat, 0);
	uint m = m_compound_alpha_size;
	for (uint k = 0; k < m_nfeat; ++k)
		{
		uint fi = m_nfeat - k - 1;
		assert(fi < m_nfeat);
		uint axis = m_axes[fi];
		byte component_code = compound_code/axis;
		codes[fi] = component_code;
		compound_code -= component_code*axis;
		}
	}

float flat_features::get_compound_subst_score_slow(
	uint8_t code1, uint8_t code2)
	{
	vector<uint8_t> code1s;
	vector<uint8_t> code2s;
	compound_code_to_component_codes(code1, code1s);
	compound_code_to_component_codes(code2, code2s);
#if DEBUG
	{
	vector<uint8_t> code1s_check;
	vector<uint8_t> code2s_check;
	uint8_t code1_check = component_codes_to_compound_code(code1s);
	uint8_t code2_check = component_codes_to_compound_code(code2s);
	assert(code1_check == code1);
	assert(code2_check == code2);
	}
#endif
	float score = 0;
	for (uint32_t fi = 0; fi < m_nfeat; ++fi)
		{
		const uint32_t AS_fi = m_alpha_sizes[fi];
		const float *logodds_fi = m_weighted_logoddsvec[fi];
		const float *logodds_row = logodds_fi + code1s[fi]*AS_fi;
		score += logodds_row[code2s[fi]];
		}
	return score;
	}

void flat_features::get_compound_logodds_slow(vector<float> &logodds)
	{
	uint compound_alpha_size = get_compound_alpha_size();
	logodds.clear();
	logodds.resize(compound_alpha_size*compound_alpha_size, FLT_MAX);
	const uint nfeat = get_nfeat();
	for (uint compound_code1 = 0; compound_code1 < compound_alpha_size;
		++compound_code1)
		{
		for (uint compound_code2 = 0; compound_code2 < compound_alpha_size;
			++compound_code2)
			{
			float score = get_compound_subst_score_slow(
				compound_code1, compound_code2);
			logodds[compound_code1*compound_alpha_size + compound_code2] = score;
			}
		}

// check symmetry
	for (uint compound_code1 = 0; compound_code1 < compound_alpha_size;
		++compound_code1)
		{
		for (uint compound_code2 = 0; compound_code2 < compound_alpha_size;
			++compound_code2)
			{
			float score12 = logodds[compound_code1*compound_alpha_size + compound_code2];
			float score21 = logodds[compound_code2*compound_alpha_size + compound_code1];
			asserta(feq(score12, score21));
			}
		}
	}

void flat_features::load_alphas(
	const vector<string> &feature_names,
	const string &logoddsfnpattern)
	{

	flat_features::init(feature_names);
	flat_features::read_logoddsvec_pattern(logoddsfnpattern);
	flat_features::set_feature_block_offsets();
	flat_features::set_symbolsvec();
	}
