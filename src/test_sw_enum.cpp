#include "myutils.h"
#include "xdpmem.h"
#include "getticks.h"
#include "flat_helpers.h"

static uint32_t s_nfeat = 3;
static uint32_t s_minL = 3;
static uint32_t s_maxL = 10;
static uint32_t s_nprof = 10;

static float s_open = -1;
static float s_ext = -0.1f;

float SWFast_Callback(XDPMem &Mem,
	uint LA,
	uint LB,
	colscorefn sf,
	float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
	string &Path);

static uint8_t **s_profs;
static uint32_t *s_prof_lengths;

static uint32_t *s_alpha_sizes;
static uint32_t s_sum_alpha_sizes;
static uint32_t *s_feature_block_offsets;
static float **s_weighted_logoddsmxvec;

static uint8_t *make_random_profile(uint32_t L)
	{
	uint8_t *prof = myalloc(uint8_t, s_nfeat*L);
	uint32_t k = 0;
	for (uint fi = 0; fi < s_nfeat; ++fi)
		{
		uint32_t alpha_size = s_alpha_sizes[fi];
		for (uint pos = 0; pos < L; ++pos)
			prof[k++] = uint8_t(randu32()%alpha_size);
		}
	assert(k == s_nfeat*L);
	return prof;
	}

static float get_random_score()
	{
	int i = int(randu32()%10) - 3;
	float r = float(randu32()%1000 + 1)/3000.0f;
	return float(i) + r;
	}

static float get_random_positive_score()
	{
	uint i = randu32()%10 + 1;
	assert(i > 0);
	float r = float(randu32()%1000 + 1)/3000.0f;
	return float(i) + r;
	}

static float *make_random_logoddsmx(uint32_t alpha_size)
	{
	float *mx = myalloc(float, alpha_size*alpha_size);
	uint n = 100;
	for (uint i = 0; i < alpha_size; ++i)
		{
		for (uint j = 0; j < alpha_size; ++j)
			{
			float score = (i == j) ? 
				get_random_positive_score() :
				get_random_score();
			mx[i*alpha_size + j] = score;
			}
		}
	return mx;
	}

static uint s_i;
static uint s_j;

static uint32_t s_L_i;
static uint32_t s_L_j;

static uint8_t *s_prof_i;
static uint8_t *s_prof_j;

static float *s_pssm_j;

static void cache_prof_j(uint j)
	{
	s_j = j;
	s_L_j = s_prof_lengths[j];
	s_prof_j = s_profs[j];
	fill_flat_pssm(s_prof_j, s_L_j, s_nfeat, s_alpha_sizes,
		s_feature_block_offsets, s_weighted_logoddsmxvec, s_pssm_j);
	}

static float prof_col_score(uint pos_i, uint pos_j)
	{
	assert(pos_i < s_L_i);
	assert(pos_j < s_L_j);
	float score = 0;
	for (uint32_t fi = 0; fi < s_nfeat; ++fi)
		{
		const uint32_t AS_fi = s_alpha_sizes[fi];
		const float *weighted_logoddsmx_fi = s_weighted_logoddsmxvec[fi];
		const uint8_t *prof_i_fi = s_prof_i + fi*s_L_i;
		const uint8_t *prof_j_fi = s_prof_j + fi*s_L_j;
		const uint8_t code_i_pos_i = prof_i_fi[pos_i];
		const float *weighted_logoddsmx_row = weighted_logoddsmx_fi + code_i_pos_i*AS_fi;
		const uint8_t code_j_pos_j = prof_j_fi[pos_j];
		score += weighted_logoddsmx_row[code_j_pos_j];
		}
	return score;
	}

static float score_path(uint start_i, uint start_j,
	const char *path, uint ncol)
	{
	uint pos_i = start_i;
	uint pos_j = start_j;
	float score = 0;
	for (uint col = 0; col < ncol; ++col)
		{
		switch (path[col])
			{
		case 'M':
			score += prof_col_score(pos_i, pos_j);
			++pos_i;
			++pos_j;
			break;
		
		case 'D':
			assert(col > 0);
			score += (path[col-1] == 'M') ? s_open : s_ext;
			++pos_i;
			break;

		case 'I':
			assert(col > 0);
			score += (path[col-1] == 'M') ? s_open : s_ext;
			++pos_j;
			break;

		default:
			Die("path[%u]='%c'", col, path[col]);
			}
		}
	assert(pos_i <= s_L_i);
	assert(pos_j <= s_L_j);
	return score;
	}

static void align_prof_enum()
	{
	void enum_sw_paths(
		uint32_t LA,
		uint32_t LB,
		vector<uint32_t>& starts_A,
		vector<uint32_t>& starts_B,
		vector<string>& paths);

	vector<uint32_t> starts_i;
	vector<uint32_t> starts_j;
	vector<string> paths;
	enum_sw_paths(s_L_i, s_L_j, starts_i, starts_j, paths);
	const uint npath = SIZE(paths);
	float best_score = 0;
	uint best_start_i = 0;
	uint best_start_j = 0;
	string best_path;
	for (uint pathidx = 0; pathidx < npath; ++pathidx)
		{
		uint start_i = starts_i[pathidx];
		uint start_j = starts_j[pathidx];
		const string &path = paths[pathidx];
		const char *path_buffer = path.c_str();
		uint ncol = uint(path.size());
		float score = score_path(start_i, start_j,
			path_buffer, ncol);
		if (score > best_score)
			{
			best_score = score;
			best_start_i = start_i;
			best_start_j = start_j;
			best_path = path;
			}
		}
	Log("%10.3g  %s  enum\n", best_score, best_path.c_str());
	}

static void align_prof_callback()
	{
	XDPMem Mem;
	uint Loi, Loj, Leni, Lenj;
	string Path;
	float score = SWFast_Callback(Mem, s_L_i, s_L_j, prof_col_score,
		s_open, s_ext,
		Loi, Loj, Leni, Lenj,
		Path);
	const char *path_buffer = Path.c_str();
	uint ncol = uint(Path.size());
	float score2 = score_path(Loi, Loj, path_buffer, ncol);
	Log("%10.3g  %s  callback\n", score, path_buffer);
	if (!feq(score, score2))
		Die("callback %.3g %.3g", score, score2);
	}

static float *s_scratch_rows;
static uint8_t *s_TB;
static void align_prof_flat()
	{
	XDPMem Mem;
	uint Loi, Loj;
	uint ncol;
	char *path_buffer = myalloc(char, 2*s_maxL);
	float score = sw_flat(s_scratch_rows, s_TB, s_L_i, s_L_j,
		prof_col_score, s_open, s_ext,
		Loi, Loj, path_buffer, ncol);
	float score2 = score_path(Loi, Loj, path_buffer, ncol);
	Log("%10.3g  %s  flat\n", score, path_buffer);
	if (!feq(score, score2))
		Die("flat %.3g %.3g", score, score2);
	myfree(path_buffer);
	}

static const float **__restrict s_scratch_pssms;
static void align_prof_pssm()
	{
	uint Loi, Loj;
	uint ncol;
	char *path_buffer = myalloc(char, 2*s_maxL);
	float score = sw_flat_pssm(
		s_scratch_rows, s_TB, s_scratch_pssms,
		s_prof_i, s_L_i,
		s_pssm_j, s_L_j, s_feature_block_offsets,
		s_nfeat, s_open, s_ext,
		Loi, Loj, path_buffer, ncol);
	Log("%10.3g  %s  pssm\n", score, path_buffer);
	float score2 = score_path(Loi, Loj, path_buffer, ncol);
	if (!feq(score, score2))
		Die("pssm %.3g %.3g", score, score2);
	myfree(path_buffer);
	}

static void align_prof_i(uint i)
	{
	s_i = i;
	s_L_i = s_prof_lengths[i];
	s_prof_i = s_profs[i];

	Log("\n======= (%u, %u)\n", s_i, s_j);
	align_prof_callback();
	align_prof_flat();
	align_prof_enum();
	align_prof_pssm();
	}

void cmd_test_sw_enum()
	{
	s_alpha_sizes = myalloc(uint32_t, s_nfeat);
	s_alpha_sizes[0] = 3;
	s_alpha_sizes[1] = 4;
	s_alpha_sizes[2] = 5;

	s_feature_block_offsets = myalloc(uint32_t, s_nfeat);
	s_scratch_rows = myalloc(float, 2*s_maxL + 2);
	s_scratch_pssms = myalloc(const float *, s_nfeat);
	s_TB = myalloc(uint8_t, s_maxL*s_maxL);

	s_sum_alpha_sizes = get_flat_pssm_feature_block_offsets(
		s_nfeat, s_alpha_sizes, s_feature_block_offsets);
	s_pssm_j = myalloc(float, s_maxL * s_sum_alpha_sizes);

	s_weighted_logoddsmxvec = myalloc(float *, s_nfeat);
	for (uint fi = 0; fi < s_nfeat; ++fi)
		{
		uint alpha_size = s_alpha_sizes[fi];
		s_weighted_logoddsmxvec[fi] = make_random_logoddsmx(alpha_size);
		}

	s_profs = myalloc(uint8_t *, s_nprof);
	s_prof_lengths = myalloc(uint32_t, s_nprof);
	for (uint i = 0; i < s_nprof; ++i)
		{
		uint32_t L = s_minL + randu32()%(s_maxL - s_minL);
		s_prof_lengths[i] = L;
		s_profs[i] = make_random_profile(L);
		}

	for (uint j = 0; j < s_nprof; ++j)
		{
		cache_prof_j(j);
		for (uint i = 0; i < s_nprof; ++i)
			align_prof_i(i);
		}
	}
