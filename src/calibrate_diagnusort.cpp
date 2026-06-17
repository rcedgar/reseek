#include "myutils.h"
#include "lookup.h"
#include "fastbench.h"
#include "triangle.h"
#include <set>

static FastBench s_FB_storage;
static FastBench *s_FB = 0;

static uint s_nhit = 0;

static vector<uint> s_query_domidx;
static vector<uint> s_target_domidx;
static vector<bool> s_is_tp;

static vector<float> s_TS;
static vector<float> s_nu_fwd_score;
static vector<float> s_nu_combined_score;
static vector<uint> s_kappa_diag_score;

static vector<vector<uint> > s_target_domidx_to_hitidxs;

static void read_hits_tsv(const string &fn, const lookup &look)
	{
	s_nhit = 0;
	s_query_domidx.clear();
	s_target_domidx.clear();
	s_is_tp.clear();
	s_TS.clear();
	s_nu_fwd_score.clear();
	s_nu_combined_score.clear();
	s_kappa_diag_score.clear();

	set<string> skipped_domains;
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
		asserta(flds.size() >= 6);

		const string &q = flds[0];
		const string &t = flds[1];
		uint qdomidx = look.get_domidx(q, true);
		uint tdomidx = look.get_domidx(t, true);
		if (qdomidx == UINT_MAX)
			skipped_domains.insert(q);
		if (tdomidx == UINT_MAX)
			skipped_domains.insert(t);
		if (qdomidx == UINT_MAX || tdomidx == UINT_MAX)
			continue;
		if (qdomidx == tdomidx)
			continue;
		if (look.is_ignored_ij(qdomidx, tdomidx))
			continue;

		s_query_domidx.push_back(qdomidx);
		s_target_domidx.push_back(tdomidx);
		s_is_tp.push_back(look.is_tp_ij(qdomidx, tdomidx));
		s_TS.push_back(StrToFloatf(flds[2]));
		s_nu_fwd_score.push_back(StrToFloatf(flds[3]));
		s_nu_combined_score.push_back(StrToFloatf(flds[4]));
		s_kappa_diag_score.push_back(StrToUint(flds[5]));
		++nline;
		}
	CloseStdioFile(f);
	ProgressFileDone();

	if (nline == 0)
		Die("No valid rows in %s", fn.c_str());
	if (skipped_domains.size() > 0)
		ProgressLog("%u domains not in lookup\n",
			uint(skipped_domains.size()));

	s_nhit = nline;

	s_target_domidx_to_hitidxs.clear();
	s_target_domidx_to_hitidxs.resize(look.get_ndom());
	for (uint hitidx = 0; hitidx < s_nhit; ++hitidx)
		{
		uint tdomidx = s_target_domidx[hitidx];
		s_target_domidx_to_hitidxs[tdomidx].push_back(hitidx);
		}

	ProgressLog("%u hits, %u domains in lookup\n",
		s_nhit, look.get_ndom());
	}

static void fill_scores_from_TS()
	{
	asserta(s_FB);
	asserta(s_FB->m_Scores);
	const uint npair = s_FB->m_npair;
	const bool top_mode =
		(s_FB->m_look->m_LT == LT_TOP_SF ||
		 s_FB->m_look->m_LT == LT_TOP_FOLD);
	const float absent =
		top_mode ? FLT_MAX : s_FB->get_missing_score();
	for (uint k = 0; k < npair; ++k)
		s_FB->m_Scores[k] = absent;

	const uint ndom = s_FB->m_look->get_ndom();
	for (uint hitidx = 0; hitidx < s_nhit; ++hitidx)
		{
		uint pairk = triangle_ij_to_k2(
			s_query_domidx[hitidx],
			s_target_domidx[hitidx],
			ndom);
		asserta(pairk < npair);
		s_FB->m_Scores[pairk] = s_TS[hitidx];
		}
	}

static double Eval3(const vector<string> &xv)
	{
	asserta(s_FB);
	(void)xv;

	s_FB->ClearHitsAndResults();
	fill_scores_from_TS();

	if (s_FB->m_look->m_LT != LT_TOP_SF &&
	    s_FB->m_look->m_LT != LT_TOP_FOLD)
		s_FB->SetScoreOrder();

	return s_FB->Bench();
	}

void cmd_calibrate_diagnusort()
	{
	const string &hits_fn = g_Arg1;
	const string lookup_fn =
		optset_lookup ? opt(lookup) : "../data/scop40x.lookup";

	s_FB = &s_FB_storage;
	s_FB->m_scores_are_evalues = opt(scores_are_evalues);
	s_FB->ReadLookup(lookup_fn);
	if (optset_dope)
		s_FB->ReadDope(opt(dope));
	s_FB->Alloc();

	read_hits_tsv(hits_fn, *s_FB->m_look);

	vector<string> xv;
	double Sum3 = Eval3(xv);
	ProgressLog("Eval3 Sum3=%.3f\n", Sum3);
	}
