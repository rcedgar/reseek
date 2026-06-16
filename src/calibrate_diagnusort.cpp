#include "myutils.h"
#include "lookup.h"
#include <set>

static lookup s_look;
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

void cmd_calibrate_diagnusort()
	{
	const string &hits_fn = g_Arg1;
	const string lookup_fn =
		optset_lookup ? opt(lookup) : "../data/scop40x.lookup";

	s_look.from_tsv(lookup_fn);
	read_hits_tsv(hits_fn, s_look);
	}
