#include "myutils.h"
#include "top_bench.h"
#include "sort.h"

void trunc_label(const string &Label,
	string &TruncatedLabel);

void top_bench::read_lookup(const string &arg_fn)
	{
	if (!optset_truth)
		{
		optset_truth = true;
		opt_truth = mystrsave("sfx");
		}
	asserta(m_look == 0);
	const string &fn =
		(arg_fn == "" ? "../data/scop40x.lookup" : arg_fn);
	m_scores_are_evalues = opt(scores_are_evalues);
	m_look = new lookup;
	m_look->from_tsv(fn);
	}

void top_bench::alloc()
	{
	if (m_score_top_tp != 0) return;
	asserta(m_look != 0);
	const uint ndom = m_look->get_ndom();
	m_score_top_tp = myalloc(float, ndom);
	m_score_top_fp = myalloc(float, ndom);
	m_domidx_top_tp = myalloc(uint, ndom);
	m_domidx_top_fp = myalloc(uint, ndom);
	m_scores = myalloc(float, 2*ndom);
	m_tps = myalloc(bool, 2*ndom);
	m_qs = myalloc(uint, 2*ndom);
	m_ts = myalloc(uint, 2*ndom);
	m_order = myalloc(uint, 2*ndom);
	}

void top_bench::clear_hits_and_results()
	{
	asserta(m_score_top_tp != 0);
	asserta(m_score_top_fp != 0);
	const uint ndom = m_look->get_ndom();
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		m_score_top_tp[domidx] = FLT_MAX;
		m_score_top_fp[domidx] = FLT_MAX;
		m_domidx_top_tp[domidx] = UINT_MAX;
		m_domidx_top_fp[domidx] = UINT_MAX;
		}
	m_topsum3 = FLT_MAX;
	m_S0_001 = FLT_MAX;
	m_S0_01 = FLT_MAX;
	m_S0_1 = FLT_MAX;
	}

void top_bench::read_tophits(const string &fn)
	{
	Die("TODO");
	FILE *f = OpenStdioFile(fn);
	string line;
	vector<string> flds;
	const uint ndom = m_look->get_ndom();
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		bool ok = ReadLineStdioFile(f, line);
		asserta(ok);
		Split(line, flds, '\t');
		asserta(flds.size() == 5);
		const string &qlabel = flds[0];
		const string &tdom_tp = flds[1];
		const string &sscore_tp = flds[2];
		const string &tdom_fp = flds[3];
		const string &sscore_fp = flds[4];

		string qdom;
		trunc_label(qlabel, qdom);
		asserta(qdom == m_look->m_doms[domidx]);

		if (sscore_tp == ".")
			{
			m_score_top_tp[domidx] = FLT_MAX;
			m_domidx_top_tp[domidx] = UINT_MAX;
			}
		else
			{
			m_score_top_tp[domidx] = StrToFloatf(sscore_tp);
			m_domidx_top_tp[domidx] = m_look->get_domidx(tdom_tp);
			}

		if (sscore_fp == ".")
			{
			m_score_top_fp[domidx] = FLT_MAX;
			m_domidx_top_fp[domidx] = UINT_MAX;
			}
		else
			{
			m_score_top_fp[domidx] = StrToFloatf(sscore_fp);
			m_domidx_top_fp[domidx] = m_look->get_domidx(tdom_fp);
			}
		}
	CloseStdioFile(f);
	}

void top_bench::read_hits(
	const string &FN,
	uint qidx,
	uint tidx,
	uint scoreidx,
	bool triangle)
	{
	GetStemName(FN, m_name);
	const uint maxidx = max(max(qidx, tidx), scoreidx);

	alloc();
	clear_hits_and_results();
	const uint ndom = m_look->get_ndom();

	FILE *f = OpenStdioFile(FN);
	string line;
	vector<string> flds;
	uint ntp = 0;
	uint n = 0;
	uint counter = 0;
	uint ntp_dope = 0;
	uint nfp_dope = 0;
	uint64 FileSize = GetStdioFileSize64(f);
	time_t lastt = time(0);
	set<string> missing;
	while (ReadLineStdioFile(f, line))
		{
		if (StartsWith("line", "# ") || StartsWith(line, "q"))
			continue;
		if (++counter%100000 == 0)
			{
			time_t t = time(0);
			if (t - lastt > 0)
				{
				uint64 FilePos = GetStdioFilePos64(f);
				double Pct = FilePos*100.0/FileSize;
				Progress("Hits %.1f%%\r", Pct);
				lastt = t;
				}
			}
		Split(line, flds, '\t');
		asserta(flds.size() > maxidx);
		const string &q = flds[qidx];
		if (q == "query") continue;
		const string &t = flds[tidx];
		uint domidxq = m_look->get_domidx(q, true);
		uint domidxt = m_look->get_domidx(t, true);
		if (domidxq == UINT_MAX)
			missing.insert(q);
		if (domidxt == UINT_MAX)
			missing.insert(t);
		if (domidxq == UINT_MAX || domidxt == UINT_MAX)
			continue;
		if (domidxq == domidxt)
			continue;
		bool ignore = is_ignored(domidxq, domidxt);
		if (ignore) continue;

		float score = (float) StrToFloat(flds[scoreidx]);

		const float top_score_tpq = m_score_top_tp[domidxq];
		const float top_score_fpq = m_score_top_fp[domidxq];
		bool tp = is_tp(domidxq, domidxt);
		bool fp = is_fp(domidxq, domidxt);
		asserta(!(tp && fp));

		if (tp)
			{
			if (top_score_tpq == FLT_MAX || better(score, top_score_tpq))
				{
				m_score_top_tp[domidxq] = score;
				m_domidx_top_tp[domidxq] = domidxt;
				}
			}
		if (fp)
			{
			if (top_score_fpq == FLT_MAX || better(score, top_score_fpq))
				{
				m_score_top_fp[domidxq] = score;
				m_domidx_top_fp[domidxq] = domidxt;
				}
			}
		if (!triangle) continue;

		const float top_score_tpt = m_score_top_tp[domidxt];
		const float top_score_fpt = m_score_top_fp[domidxt];
		if (tp)
			{
			if (top_score_tpt == FLT_MAX || better(score, top_score_tpt))
				{
				m_score_top_tp[domidxt] = score;
				m_domidx_top_tp[domidxt] = domidxq;
				}
			}
		else
			{
			if (top_score_fpt == FLT_MAX || better(score, top_score_fpt))
				{
				m_score_top_fp[domidxt] = score;
				m_domidx_top_fp[domidxt] = domidxq;
				}
			}
		}

	Progress("Hits 100.00%%\n");
	const uint nmiss = uint(missing.size());
	if (nmiss > 0)
		{
		ProgressLog("%u domains not in lookup\n", uint(nmiss));
		uint n = 0;
		for (auto dom : missing)
			{
			Log(">%s\n", dom.c_str());
			if (++n > 10)
				{
				Log("... %u more\n", nmiss - n);
				break;
				}
			}
		}
	CloseStdioFile(f);
	}

double top_bench::bench(const string &fn, const string &msg)
	{
	FILE *f = 0;
	if (fn != "")
		f = CreateStdioFile(fn);

	asserta(m_score_top_tp != 0);
	asserta(m_score_top_fp != 0);
	const uint ndom = m_look->get_ndom();
	const uint nsingle = m_look->get_singleton_count();
	const uint non_singleton_count = ndom - nsingle;
	
	m_S0_001 = FLT_MAX;
	m_S0_01 = FLT_MAX;
	m_S0_1 = FLT_MAX;

	uint n = 0;
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		float score_tp = m_score_top_tp[domidx];
		float score_fp = m_score_top_fp[domidx];
		if (score_tp != FLT_MAX)
			{
			m_scores[n] = score_tp;
			m_tps[n] = true;
			m_qs[n] = domidx;
			m_ts[n] = m_domidx_top_tp[domidx];
			++n;
			}

		if (score_fp != FLT_MAX)
			{
			m_scores[n] = score_fp;
			m_tps[n] = false;
			m_qs[n] = domidx;
			m_ts[n] = m_domidx_top_fp[domidx];
			++n;
			}
		}
	const uint N = n;

	if (m_scores_are_evalues)
		QuickSortOrder<float>(m_scores, N, m_order);
	else
		QuickSortOrderDesc<float>(m_scores, N, m_order);

	uint nt_top = 0;
	uint nf_top = 0;
	float top_EPQ = 0;
	float top_Sens = 0;
	float last_score = m_scores[m_order[0]];
	for (uint orderidx = 0; orderidx < N; ++orderidx)
		{
		uint scoreidx = m_order[orderidx];
		float score = m_scores[scoreidx];
		asserta(score != FLT_MAX);
		bool tp = m_tps[scoreidx];
		if (score != last_score)
			{
			asserta(better(last_score, score));
			if (top_EPQ <= 0.001) m_S0_001 = top_Sens;
			if (top_EPQ <= 0.01) m_S0_01   = top_Sens;
			if (top_EPQ <= 0.1)  m_S0_1    = top_Sens;

			last_score = score;
			}
		if (tp)
			++nt_top;
		else
			++nf_top;

		top_EPQ = float(nf_top)/ndom;
		top_Sens = float(nt_top)/non_singleton_count;

		if (f == 0) continue;
		uint domidxq = m_qs[scoreidx];
		uint domidxt = m_ts[scoreidx];
		uint famidxq = m_look->m_domidx2famidx[domidxq];
		uint famidxt = m_look->m_domidx2famidx[domidxt];
		const string &q = m_look->get_dom(domidxq);
		const string &t = m_look->get_dom(domidxt);
		fprintf(f, "%s/%s", q.c_str(), m_look->m_fams[famidxq].c_str());
		fprintf(f, "\t%s/%s", t.c_str(), m_look->m_fams[famidxt].c_str());
		fprintf(f, "\t%.3g", m_scores[scoreidx]);
		fprintf(f, "\t%c", m_tps[scoreidx] ? 'T' : 'F');
		fprintf(f, "\t%.3g", top_Sens);
		fprintf(f, "\t%.3g", top_EPQ);
		fprintf(f, "\n");
		}

	m_topsum3 = m_S0_001*2 + m_S0_01*3/2 + m_S0_1;

	if (msg != "noshow")
		{
		if (msg != "")
			ProgressLog("%s ", msg.c_str());
		ProgressLog("S_0.001=%.3f", m_S0_001);
		ProgressLog(" S_0.01=%.3f", m_S0_01);
		ProgressLog(" S_0.1=%.3f", m_S0_1);
		ProgressLog(" Top3=%.3f", m_topsum3);
		ProgressLog(" %s", m_look->get_truthstr());
		if (m_name != "")
			ProgressLog(" %s", m_name.c_str());
		ProgressLog("\n");
		}

	CloseStdioFile(f);
	return m_topsum3;
	}

void cmd_top_bench_tophits()
	{
	const string &hitsfn = g_Arg1;
	const string lookupfn =
		(optset_lookup ? opt(lookup) : "../data/scop40x.lookup");

	top_bench TB;
	TB.read_lookup(lookupfn);
	TB.alloc();
	TB.read_tophits(hitsfn);
	TB.bench(opt(output));
	}

void cmd_top_bench_hits()
	{
	const string &hitsfn = g_Arg1;
	const string lookupfn =
		(optset_lookup ? opt(lookup) : "../data/scop40x.lookup");

	top_bench TB;
	TB.read_lookup(lookupfn);
	uint qidx = UINT_MAX;
	uint tidx = UINT_MAX;
	uint scoreidx = UINT_MAX;
	if (optset_qfield)
		{
		asserta(optset_tfield);
		asserta(optset_scorefield);
		qidx = opt(qfield) - 1;
		tidx = opt(tfield) - 1;
		scoreidx = opt(scorefield) - 1;
		}
	else
		guess_fields(
			hitsfn, qidx, tidx, scoreidx, TB.m_scores_are_evalues);

	ProgressLog("hits %u,%u,%u %s\n",
		qidx+1, tidx+1, scoreidx+1,
		TB.m_scores_are_evalues ? "E-values" : "scores");

	TB.read_hits(hitsfn, qidx, tidx, scoreidx, opt(triangle));
	TB.bench(opt(output));
	}

void cmd_top_bench()
	{
	}
