#if 0
#include "myutils.h"
#include "lookup.h"
#include "fastbench.h"
#include "triangle.h"
#include <algorithm>
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
static vector<vector<uint> > s_domidx_to_hitidxs_in_decreasing_diag_score_order;

// minfwd=20;nfminfwd=140;nfmincmb=70; SF
static float s_min_nu_fwd_score = 140;
static float s_min_nu_combined_score = 70;
static uint s_nparaln_fwd = 0;
static uint s_nparaln_rev = 0;
static uint s_nmegaln = 0;

static bool stop_nu()
	{
	return false;
	}

static bool stop_mega()
	{
	return false;
	}

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
	s_domidx_to_hitidxs_in_decreasing_diag_score_order.clear();

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
	s_domidx_to_hitidxs_in_decreasing_diag_score_order.resize(look.get_ndom());
	for (uint hitidx = 0; hitidx < s_nhit; ++hitidx)
		{
		uint qdomidx = s_query_domidx[hitidx];
		uint tdomidx = s_target_domidx[hitidx];
		s_target_domidx_to_hitidxs[tdomidx].push_back(hitidx);
		s_domidx_to_hitidxs_in_decreasing_diag_score_order[qdomidx].push_back(hitidx);
		}

	for (uint qdomidx = 0; qdomidx < look.get_ndom(); ++qdomidx)
		{
		vector<uint> &hitidxs = s_domidx_to_hitidxs_in_decreasing_diag_score_order[qdomidx];
		if (hitidxs.empty())
			continue;
		std::sort(hitidxs.begin(), hitidxs.end(),
			[](uint a, uint b)
				{
				return s_kappa_diag_score[a] > s_kappa_diag_score[b];
				});
		}

	ProgressLog("%u hits, %u domains in lookup\n",
		s_nhit, look.get_ndom());
	}

static void set_all_scores_to_missing()
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
	}

static void fill_scores_from_TS()
	{
	set_all_scores_to_missing();

	const uint ndom = s_FB->m_look->get_ndom();
	const uint npair = s_FB->m_npair;
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

// fills all scores in diag order, the ordering has
// no effect on the result, this is a skeleton
// for builing functional variants.
static void fill_scores_diag_order()
	{
	set_all_scores_to_missing();

	const uint ndom = s_FB->m_look->get_ndom();
	const uint npair = s_FB->m_npair;
	for (uint qdomidx = 0; qdomidx < ndom; ++qdomidx)
		{
		const vector<uint> &hitidxs =
			s_domidx_to_hitidxs_in_decreasing_diag_score_order[qdomidx];
		for (uint i = 0; i < uint(hitidxs.size()); ++i)
			{
			const uint hitidx = hitidxs[i];
			uint pairk = triangle_ij_to_k2(
				s_query_domidx[hitidx],
				s_target_domidx[hitidx],
				ndom);
			asserta(pairk < npair);
			s_FB->m_Scores[pairk] = s_TS[hitidx];
			}
		}
	}

static void fill_scores(bool do_nu_stop, bool do_mega_stop)
	{
	set_all_scores_to_missing();

	const uint ndom = s_FB->m_look->get_ndom();
	const uint npair = s_FB->m_npair;
	s_nparaln_fwd = 0;
	s_nparaln_rev = 0;
	s_nmegaln = 0;
	for (uint qdomidx = 0; qdomidx < ndom; ++qdomidx)
		{
		const vector<uint> &hitidxs =
			s_domidx_to_hitidxs_in_decreasing_diag_score_order[qdomidx];
		for (uint i = 0; i < uint(hitidxs.size()); ++i)
			{
			const uint hitidx = hitidxs[i];
			asserta(hitidx < s_kappa_diag_score.size());
			uint diag_score = s_kappa_diag_score[hitidx];

			++s_nparaln_fwd;
			float nu_fwd_score = s_nu_fwd_score[hitidx];
			if (nu_fwd_score < s_min_nu_fwd_score) continue;

			++s_nparaln_rev;
			float nu_combined_score = s_nu_combined_score[hitidx];
			if (nu_combined_score < s_min_nu_combined_score) continue;

			if (do_nu_stop && stop_nu())
				break;

			uint pairk = triangle_ij_to_k2(
				s_query_domidx[hitidx],
				s_target_domidx[hitidx],
				ndom);
			asserta(pairk < npair);
			++s_nmegaln;
			if (do_mega_stop && stop_mega())
				break;

			s_FB->m_Scores[pairk] = s_TS[hitidx];
			}
		}
	}

static double Eval3_all_TS(const vector<string> &xv)
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

static double Eval3_diag_order(const vector<string> &xv)
	{
	asserta(s_FB);
	(void)xv;

	s_FB->ClearHitsAndResults();
	fill_scores_diag_order();

	if (s_FB->m_look->m_LT != LT_TOP_SF &&
	    s_FB->m_look->m_LT != LT_TOP_FOLD)
		s_FB->SetScoreOrder();

	return s_FB->Bench();
	}

static double Eval3(const vector<string> &xv,
	bool do_diag_stop, bool do_nu_stop)
	{
	asserta(s_FB);
	(void)xv;

	s_FB->ClearHitsAndResults();
	fill_scores(do_diag_stop, do_nu_stop);

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
	uint nhit = uint(s_TS.size());

	vector<string> xv;
	double Sum3_0 = Eval3(xv, false, false);
	double t_0 = s_nparaln_fwd + s_nparaln_rev + 50*s_nmegaln;
	ProgressLog("Eval3(false, false) Sum3=%.3f nfwd=%u nrev=%u nmega=%u t=%.3g\n",
		Sum3_0, s_nparaln_fwd, s_nparaln_rev, s_nmegaln, t_0);

	double Sum3 = Eval3(xv, true, true);
	double t = s_nparaln_fwd + s_nparaln_rev + 50*s_nmegaln;

	ProgressLog("Eval3(true, true) Sum3=%.3f nfwd=%u nrev=%u nmega=%u t=%.3g\n",
		Sum3, s_nparaln_fwd, s_nparaln_rev, s_nmegaln, t);
	}

#endif
