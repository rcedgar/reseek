#include "myutils.h"
#include "flat_helpers.h"
#include "flat_alphas.h"
#include "tabbedlines.h"
#include "chaq.h"
#include "fan.h"

/***
# C:\src\reseek\src\Release\reseek.exe -flat_quantize ../data/scop40x.bca ...
# [2544ec1-dirty] 2026-05-07
fan	turnd
alpha_size	32
median	1540
thresholds	31	344	520	647	750	854	...
***/

static vector<FAN> s_fans;
static vector<uint> s_alpha_sizes;
static vector<uint16_t> s_medians;
static vector<uint16_t *> s_thresholds;

void flat_alphas::init_from_alphadir(
	const string &arg_alphadir,
	const vector<string> &alpha_names)
	{
	asserta(!alpha_names.empty());

	string alphadir = arg_alphadir;
	Dirize(alphadir);

	s_fans.clear();
	s_alpha_sizes.clear();
	s_medians.clear();
	s_thresholds.clear();

	const uint nfeat = uint(alpha_names.size());
	alloc(nfeat);
	m_alpha_names = alpha_names;
	m_sum_alpha_sizes = 0;
	m_compound_alpha_size = 1;
	m_entropyfi = UINT_MAX;

	string compound;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		const string &alpha_name = alpha_names[fi];
		if (fi > 0)
			compound += "+";
		compound += alpha_name;
		uint alpha_size = 0;
		FAN fan = parse_alpha_name(alpha_name, alpha_size);
		m_fans.push_back(fan);

		m_alpha_sizes[fi] = alpha_size;
		if (StartsWith(alpha_name, "sec") || alpha_name == "Conf")
			m_entropyfi = fi;

		m_axes[fi] = m_compound_alpha_size;
		m_compound_alpha_size *= alpha_size;
		m_sum_alpha_sizes += alpha_size;

		string logoddsfn;
		Ps(logoddsfn, "%s/%s%u.logodds",
			alphadir.c_str(), FAN2str(fan), alpha_size);

		vector<float> logodds;
		uint alpha_size2 = read_logodds(logoddsfn, logodds);
		asserta(alpha_size2 == alpha_size);

		const uint n = alpha_size*alpha_size;
		m_unweighted_logoddsvec[fi] = myalloc(float, n);
		m_weighted_logoddsvec[fi] = myalloc(float, n);
		for (uint k = 0; k < n; ++k)
			{
			const float score = logodds[k];
			asserta(score >= MIN_SANE_SCORE && score <= MAX_SANE_SCORE);
			m_unweighted_logoddsvec[fi][k] = score;
			m_weighted_logoddsvec[fi][k] = BAD_SCORE;
			}

		uint16_t median = UINT16_MAX;
		uint16_t *thresholds = 0;
		if (is_quantized(fan))
			{
			string quantizefn;
			Ps(quantizefn, "%s/%s%u.quantize",
				alphadir.c_str(), FAN2str(fan), alpha_size);
			thresholds = read_quantize(quantizefn, alpha_size, median);
			}

		s_fans.push_back(fan);
		s_alpha_sizes.push_back(alpha_size);
		s_medians.push_back(median);
		s_thresholds.push_back(thresholds);
		}

	apply_unit_weights();
	set_feature_block_offsets();
	set_symbolsvec();

	ProgressLog("Loaded %s\n", compound.c_str());
	}

uint16_t chaq::get_undef_value(FAN fan, uint alpha_size)
	{
	const size_t n = s_fans.size();
	assert(s_alpha_sizes.size() == n);
	assert(s_medians.size() == n);
	assert(s_thresholds.size() == n);
	for (size_t i = 0; i < n; ++i)
		{
		if (s_fans[i] == fan && s_alpha_sizes[i] == alpha_size)
			{
			uint16_t median = s_medians[i];
			if (median == UINT16_MAX)
				Die("chaq::get_undef_value(%s, %u) median=UINT16_MAX",
					FAN2str(fan), alpha_size);
			return median;
			}
		}
	Die("chaq::get_undef_value(%s, %u) not found",
		FAN2str(fan), alpha_size);
	return UINT16_MAX;
	}

cp_uint16_t chaq::get_thresholds(FAN fan, uint alpha_size)
	{
	const size_t n = s_fans.size();
	assert(s_alpha_sizes.size() == n);
	assert(s_medians.size() == n);
	assert(s_thresholds.size() == n);
	for (size_t i = 0; i < n; ++i)
		{
		if (s_fans[i] == fan && s_alpha_sizes[i] == alpha_size)
			return s_thresholds[i];
		}
	Die("chaq::get_thresholds(%s, %u)", FAN2str(fan), alpha_size);
	return 0;
	}

FAN parse_alpha_name(const string &alpha_name, uint &alpha_size)
	{
	size_t n = alpha_name.size();
	asserta(n > 2);
	string ssize;
	for (int j = int(n)-1; j > 0; --j)
		{
		char c = alpha_name[j];
		if (isdigit(c))
			ssize = c + ssize;
		else
			{
			alpha_size = StrToUint(ssize);
			string fan = alpha_name.substr(0, j+1);
			return str2FAN(fan);
			}
		}
	Die("parse_alpha_name(%s)", alpha_name.c_str());
	return FAN_COUNT;
	}

uint16_t *read_quantize(const string &fn, uint alpha_size, uint16_t &median)
	{
	tabbedlines tl(fn);
	string s = tl.get_str("fan");
	FAN fan = str2FAN(s);
	uint alpha_size2 = tl.get_int("alpha_size");
	asserta(alpha_size2 == alpha_size);
	median = tl.get_int("median");
	uint16_t *thresholds = tl.get_int16_flat_vec("thresholds", alpha_size-1);
	return thresholds;
	}

static void load(const string &alphadir, FAN fan, uint alpha_size)
	{
	string logoddsfn;
	Ps(logoddsfn, "%s/%s%u.logodds",
		alphadir.c_str(), FAN2str(fan), alpha_size);

	vector<float> logodds;
	uint alpha_size2 = read_logodds(logoddsfn, logodds);
	asserta(alpha_size2 == alpha_size);

	uint16_t median = UINT16_MAX;
	uint16_t *thresholds = 0;
	if (is_quantized(fan))
		{
		string quantizefn;
		Ps(quantizefn, "%s/%s%u.quantize",
			alphadir.c_str(), FAN2str(fan), alpha_size);
		thresholds = read_quantize(quantizefn, alpha_size, median);
		}
	s_fans.push_back(fan);
	s_alpha_sizes.push_back(alpha_size);
	s_medians.push_back(median);
	s_thresholds.push_back(thresholds);
	}

void load_alphadir(const string &arg_alphadir)
	{
	string alphadir = arg_alphadir;
	Dirize(alphadir);

	vector<string> fns;
	vector<bool> subdirs;
	mylistdir(alphadir, fns, subdirs);

	const size_t n = fns.size();
	for (uint i = 0; i < n; ++i)
		{
		if (subdirs[i]) continue;
		const string &fn = fns[i];
		if (EndsWith(fn, ".logodds"))
			{
			size_t n = fn.size() - strlen(".logodds");
			string alpha_name = fn.substr(0, n);
			uint alpha_size;
			FAN fan = parse_alpha_name(alpha_name, alpha_size);
			load(alphadir, fan, alpha_size);
			}
		}
	ProgressLog("Loaded %s, %u alphabets found\n",
		alphadir.c_str(), uint(s_fans.size()));
	}

void load_alphadir_names(
	const string &arg_alphadir,
	const vector<string> &alpha_names)
	{
	string alphadir = arg_alphadir;
	Dirize(alphadir);
	const size_t n = alpha_names.size();
	for (uint i = 0; i < n; ++i)
		{
		const string &alpha_name = alpha_names[i];
		uint alpha_size;
		FAN fan = parse_alpha_name(alpha_name, alpha_size);
		load(alphadir, fan, alpha_size);
		}
	ProgressLog("Loaded %s, %u alphabets found\n",
		alphadir.c_str(), uint(s_fans.size()));
	}

void cmd_read_alphadir()
	{
	load_alphadir(g_Arg1);
	}