#include "myutils.h"
#include "fastbench.h"
#include "peaker.h"
#include <set>
#include <unordered_map>

static Peaker *s_peaker;
static FastBench *s_fb;
static float *s_data = 0;
static vector<uint> s_pair_ks;
static uint s_nfloat = 0;

static float scale_value(float v, float col_min, float col_max)
	{
	if (col_max > col_min)
		return 1 + 99*(v - col_min)/(col_max - col_min);
	return 1;
	}

static void read_feature_tsv(
	const string &fn,
	const lookup &look)
	{
	unordered_map<uint, vector<float> > k2raw;
	set<string> skipped_domains;
	uint nfloat = UINT_MAX;
	uint nline = 0;

	FILE *f = OpenStdioFile(fn);
	string line;
	vector<string> flds;
	ProgressFileInit(f, "Reading %s", fn.c_str());
	while (ReadLineStdioFile(f, line))
		{
		ProgressFileStep();
		if (line.empty())
			continue;
		Split(line, flds, '\t');
		if (nfloat == UINT_MAX)
			{
			asserta(flds.size() >= 3);
			nfloat = uint(flds.size()) - 2;
			}
		asserta(flds.size() == nfloat + 2);

		const string &q = flds[0];
		const string &t = flds[1];
		uint domidxq = look.get_domidx(q, true);
		uint domidxt = look.get_domidx(t, true);
		if (domidxq == UINT_MAX)
			skipped_domains.insert(q);
		if (domidxt == UINT_MAX)
			skipped_domains.insert(t);
		if (domidxq == UINT_MAX || domidxt == UINT_MAX)
			continue;

		vector<float> raw(nfloat);
		for (uint fi = 0; fi < nfloat; ++fi)
			{
			float v = StrToFloatf(flds[2 + fi]);
			asserta(!isnan(v));
			raw[fi] = v;
			}

		uint k = look.get_pair_idx_upper_triangle_with_diagonal(
			domidxq, domidxt);
		k2raw[k] = raw;
		++nline;
		}
	CloseStdioFile(f);
	ProgressFileDone();

	if (nfloat == UINT_MAX)
		Die("Empty feature tsv %s", fn.c_str());
	if (nline == 0)
		Die("No valid rows in %s", fn.c_str());
	if (skipped_domains.size() > 0)
		ProgressLog("%u domains not in lookup\n",
			uint(skipped_domains.size()));

	s_nfloat = nfloat;
	vector<float> col_min(nfloat, FLT_MAX);
	vector<float> col_max(nfloat, -FLT_MAX);
	for (unordered_map<uint, vector<float> >::const_iterator iter =
		k2raw.begin(); iter != k2raw.end(); ++iter)
		{
		const vector<float> &raw = iter->second;
		asserta(raw.size() == nfloat);
		for (uint fi = 0; fi < nfloat; ++fi)
			{
			float v = raw[fi];
			if (v < col_min[fi]) col_min[fi] = v;
			if (v > col_max[fi]) col_max[fi] = v;
			}
		}

	s_pair_ks.clear();
	s_pair_ks.reserve(k2raw.size());
	myfree(s_data);
	s_data = myalloc(float, uint(k2raw.size()*nfloat));

	uint pair_idx = 0;
	for (unordered_map<uint, vector<float> >::const_iterator iter =
		k2raw.begin(); iter != k2raw.end(); ++iter)
		{
		uint k = iter->first;
		const vector<float> &raw = iter->second;
		s_pair_ks.push_back(k);
		float *row = s_data + pair_idx*nfloat;
		for (uint fi = 0; fi < nfloat; ++fi)
			row[fi] = scale_value(raw[fi], col_min[fi], col_max[fi]);
		++pair_idx;
		}

	ProgressLog("%u pairs, %u features\n",
		uint(s_pair_ks.size()), nfloat);
	}

static double eval_sum3(const vector<string> &xv)
	{
	asserta(s_fb);
	asserta(s_peaker);
	vector<double> weights;
	s_peaker->xv2values(xv, weights);
	asserta(weights.size() == s_nfloat);

	const uint npair = uint(s_pair_ks.size());
	for (uint pi = 0; pi < npair; ++pi)
		{
		uint k = s_pair_ks[pi];
		const float *row = s_data + pi*s_nfloat;
		double score = 0;
		for (uint fi = 0; fi < s_nfloat; ++fi)
			score += weights[fi]*row[fi];
		asserta(!isnan(score));
		s_fb->m_Scores[k] = (float) score;
		}

	Progress("Sorting...\r");
	s_fb->SetScoreOrder_Parallel();
	return s_fb->Bench();
	}

static void optimize(
	const vector<string> &spec_lines,
	double &best_y,
	vector<string> &best_xv)
	{
	const string opt_name("latinclimb");
	string global_spec;
	Peaker::GetGlobalSpec(spec_lines, global_spec);

	uint latin_bin_count = Peaker::SpecGetInt(global_spec, "latin", UINT_MAX);
	uint hj_count = Peaker::SpecGetInt(global_spec, "hj", UINT_MAX);
	asserta(latin_bin_count != UINT_MAX);
	asserta(hj_count != UINT_MAX);

	Peaker &p = *new Peaker(0, opt_name);
	p.Init(spec_lines, eval_sum3);
	s_peaker = &p;

	ProgressLog("=========================================\n");
	ProgressLog("%s latin (%u)\n", opt_name.c_str(), latin_bin_count);
	ProgressLog("=========================================\n");
	asserta(latin_bin_count > 0);
	p.RunLatin(latin_bin_count);

	vector<uint> top_eval_idxs;
	p.GetTopEvalIdxs(hj_count, top_eval_idxs);
	const uint n = SIZE(top_eval_idxs);
	if (n == 0)
		Die("No evals %s", opt_name.c_str());

	for (uint ki = 0; ki < n; ++ki)
		{
		ProgressLog("=========================================\n");
		ProgressLog("%s HJ %u/%u\n", opt_name.c_str(), ki+1, n);
		ProgressLog("=========================================\n");

		uint eval_idx = top_eval_idxs[ki];
		string child_name;
		Ps(child_name, "HJ%u/%u", ki+1, n);
		Peaker *child = p.MakeChild(child_name);
		double y = p.m_ys[eval_idx];
		const vector<string> &xv = p.m_xvs[eval_idx];
		child->AppendResult(xv, y, "HJstart");
		child->HJ_RunHookeJeeves();
		p.AppendChildResults(*child);
		delete child;
		ProgressLog("=========================================\n");
		ProgressLog("%s HJ %u/%u converged\n", opt_name.c_str(), ki+1, n);
		ProgressLog("=========================================\n");
		}

	best_y = p.m_Best_y;
	best_xv = p.m_Best_xv;
	ProgressLog("=========================================\n");
	ProgressLog("%s completed\n", opt_name.c_str());
	ProgressLog("=========================================\n");
	}

void cmd_hjtsv()
	{
	const string &feature_fn = g_Arg1;
	const string lookup_fn =
		optset_lookup ? opt(lookup) : "../data/scop40x.lookup";

	s_fb = new FastBench;
	s_fb->ReadLookup(lookup_fn);
	s_fb->Alloc();

	read_feature_tsv(feature_fn, *s_fb->m_look);

	vector<string> spec_lines;
	spec_lines.push_back("strategy=latinclimb;");
	spec_lines.push_back("latin=32;");
	spec_lines.push_back("rates=1.3,1.05;");
	spec_lines.push_back("hj=1;");
	for (uint fi = 0; fi < s_nfloat; ++fi)
		{
		string line;
		Ps(line, "var=var%u;min=0;max=1;weight=no;sigfig=3;", fi+1);
		spec_lines.push_back(line);
		}

	double best_y;
	vector<string> best_xv;
	optimize(spec_lines, best_y, best_xv);

	string xstr;
	s_peaker->xv2xss(best_xv, xstr);
	Log("@TSV@");
	Log("\t%.4g", best_y);
	Log("\t%s", xstr.c_str());
	Log("\n");
	ProgressLog("FINAL [%.4g]  %s\n", best_y, xstr.c_str());

	vector<double> best_values;
	s_peaker->xv2values(best_xv, best_values);
	for (uint fi = 0; fi < best_values.size(); ++fi)
		ProgressLog("%s=%.4g\n",
			s_peaker->m_VarNames[fi].c_str(),
			best_values[fi]);
	}