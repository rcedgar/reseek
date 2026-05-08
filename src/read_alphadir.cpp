#include "myutils.h"
#include "flat_helpers.h"
#include "tabbedlines.h"
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
	ProgressLog("  logodds\n");

	uint16_t median = UINT16_MAX;
	uint16_t *thresholds = 0;
	if (is_quantized(fan))
		{
		string quantizefn;
		Ps(quantizefn, "%s/%s%u.quantize",
			alphadir.c_str(), FAN2str(fan), alpha_size);
		thresholds = read_quantize(quantizefn, alpha_size, median);
		ProgressLog("  quantize\n");
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
			ProgressLog("%s [%u]\n", alpha_name.c_str(), alpha_size);
			}
		}
	}

void cmd_read_alphadir()
	{
	load_alphadir(g_Arg1);
	}