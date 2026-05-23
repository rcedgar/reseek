#include "myutils.h"
#include "triangle.h"
#include "fastbench.h"
#include "sort.h"

#define SAVE_NOT_IN_DOPE	0

#if PARALLEL_SORT

#include "parallel_sort.h"

void FastBench::WriteTopHits(const string &fn) const
	{
	if (fn == "")
		return;
	FILE *f = CreateStdioFile(fn);
	const uint ndom = m_look->get_ndom();
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		const string &domq = m_look->get_dom(domidx);
		float score_top_tp = m_score_top_TP[domidx];
		float score_top_fp = m_score_top_FP[domidx];
		uint domidx_top_tp = m_domidx_top_TP[domidx];
		uint domidx_top_fp = m_domidx_top_FP[domidx];

		uint famidxq = m_look->m_domidx2famidx[domidx];

		fprintf(f, "%s/%s",
			domq.c_str(), m_look->m_fams[famidxq].c_str());

		if (domidx_top_tp == UINT_MAX)
			{
			asserta(score_top_tp == FLT_MAX);
			fprintf(f, "\t.\t.");
			}
		else
			{
			const string &domt = m_look->get_dom(domidx_top_tp);
			uint famidxt = m_look->m_domidx2famidx[domidx_top_tp];
			fprintf(f, "\t%s/%s\t%.3g",
				domt.c_str(), m_look->m_fams[famidxt].c_str(),
				score_top_tp);
			}

		if (domidx_top_fp == UINT_MAX)
			{
			asserta(score_top_fp == FLT_MAX);
			fprintf(f, "\t.\t.");
			}
		else
			{
			const string &domt = m_look->get_dom(domidx_top_fp);
			uint famidxt = m_look->m_domidx2famidx[domidx_top_fp];
			fprintf(f, "\t%s/%s\t%.3g",
				domt.c_str(), m_look->m_fams[famidxt].c_str(),
				score_top_fp);
			}
		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	}

void FastBench::SetScoreOrder_Parallel()
{
	asserta(m_Scores);
	const uint K = triangle_get_K(m_SeqCount);

	if (m_ScoreOrderCap < K)
		{
		if (m_ScoreOrder != 0)
			myfree(m_ScoreOrder);
		m_ScoreOrder = myalloc(uint, K);
		m_ScoreOrderCap = K;

		// first use: initialize identity permutation
		for (uint i = 0; i < K; ++i)
			m_ScoreOrder[i] = i;
		}

	// Reuse previous order on subsequent calls.
	// This is the important part: do NOT reinitialize every time.

	if (m_scores_are_evalues)
		QuickSortOrder_Parallel(m_Scores, K, m_ScoreOrder);
	else
		QuickSortOrderDesc_Parallel(m_Scores, K, m_ScoreOrder);
}
#endif

void FastBench::Alloc()
	{
	asserta(m_look);
	const uint ndom = m_look->get_ndom();
	const uint npair = m_look->get_pair_count_upper_triangle_with_diagonal();

	if (m_score_top_FP == 0)
		{
		asserta(m_score_top_TP == 0);
		m_score_top_TP = myalloc(float, ndom);
		m_score_top_FP = myalloc(float, ndom);
		m_domidx_top_TP = myalloc(uint, ndom);
		m_domidx_top_FP = myalloc(uint, ndom);
		}

#if PARALLEL_SORT
	if (m_Scores == 0)
		{
		m_Scores = myalloc(float, npair);
		for (uint i = 0; i < npair; ++i) m_Scores[i] = FLT_MAX;
		m_PairCount = npair;
		}
	else
		asserta(m_PairCount == npair);
#else
	m_PairCount = npair;
	myfree(m_Scores);
	myfree(m_ScoreOrder);
	m_Scores = myalloc(float, npair);
#endif
	}

void FastBench::AppendHit(uint i, uint j, float Score)
	{
	uint k = triangle_ij_to_k(i, j, m_SeqCount);
	m_Scores[k] = Score;
	SubclassAppendHit(i, j, Score);
	}

void FastBench::SetScoreOrder()
	{
#if PARALLEL_SORT
	SetScoreOrder_Parallel();
#else
	SetScoreOrder_Serial();
#endif
	}

void FastBench::SetScoreOrder_Serial()
	{
	asserta(m_Scores);
	uint K = triangle_get_K(m_SeqCount);
	if (m_ScoreOrder != 0)
		myfree(m_ScoreOrder);
	m_ScoreOrder = myalloc(uint, K);
	if (m_scores_are_evalues)
		QuickSortOrder(m_Scores, K, m_ScoreOrder);
	else
		QuickSortOrderDesc(m_Scores, K, m_ScoreOrder);
	}

bool FastBench::IsIgnored(uint LabelIdx_i, uint LabelIdx_j) const
	{
	assert(m_look);
	return m_look->is_ignored_ij(LabelIdx_i, LabelIdx_j);
	}

bool FastBench::IsTP(uint LabelIdx_i, uint LabelIdx_j) const
	{
	assert(m_look);
	return m_look->is_tp_ij(LabelIdx_i, LabelIdx_j);
	}

double FastBench::Bench(const string &Msg)
	{
	asserta(m_ScoreOrder != 0);
	const uint ndom = m_look->get_ndom();
	uint K = triangle_get_K(ndom);
	uint nt = 0;
	uint nf = 0;
	uint nt_top = 0;
	uint nf_top = 0;

	m_CVESum3 = FLT_MAX;
	m_TopSum3 = FLT_MAX;
	m_SEPQ0_1 = FLT_MAX;
	m_SEPQ1 = FLT_MAX;
	m_SEPQ10 = FLT_MAX;
	m_top_SEPQ0_001 = FLT_MAX;
	m_top_SEPQ0_01 = FLT_MAX;
	m_top_SEPQ0_1 = FLT_MAX;
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		m_score_top_TP[domidx] = FLT_MAX;
		m_score_top_FP[domidx] = FLT_MAX;
		m_domidx_top_TP[domidx] = UINT_MAX;
		m_domidx_top_FP[domidx] = UINT_MAX;
		}

	float LastScore = m_scores_are_evalues ? -9e9f : FLT_MAX;
	const uint non_singleton_count = ndom - m_look->m_singleton_count;

	m_SEPQ0_1 = FLT_MAX;
	m_SEPQ1 = FLT_MAX;
	m_SEPQ10 = FLT_MAX;

	m_top_SEPQ0_001 = FLT_MAX;
	m_top_SEPQ0_01 = FLT_MAX;

	bool triangle = opt(triangle);
	for (uint k = 0; k < K; ++k)
		{
		uint HitIdx = m_ScoreOrder[k];
		uint domidx_i, domidx_j;
		triangle_k_to_ij(HitIdx, ndom, domidx_i, domidx_j);
		if (domidx_i == domidx_j)
			continue;
		float Score = m_Scores[HitIdx];
		if (Score == FLT_MAX) continue;
		if (Score != LastScore)
			{
			if (m_scores_are_evalues)
				asserta(Score > LastScore);
			else
				asserta(Score < LastScore);
			float EPQ = 2*float(nf)/ndom;
			float Sens = 2*float(nt)/m_look->m_NT;
			if (m_SEPQ0_1 == FLT_MAX && EPQ >= 0.1) m_SEPQ0_1 = Sens;
			if (m_SEPQ1 == FLT_MAX   && EPQ >= 1)   m_SEPQ1   = Sens;
			if (m_SEPQ10 == FLT_MAX  && EPQ >= 10)  m_SEPQ10  = Sens;

			LastScore = Score;
			}
		if (IsIgnored(domidx_i, domidx_j))
			continue;

		if (IsTP(domidx_i, domidx_j))
			{
			++nt;
			if (m_score_top_TP[domidx_i] == FLT_MAX ||
				Score > m_score_top_TP[domidx_i])
				{
				m_score_top_TP[domidx_i] = Score;
				m_domidx_top_TP[domidx_i] = domidx_j;
				}
			if (triangle)
				{
				if (m_score_top_TP[domidx_j] == FLT_MAX ||
					Score > m_score_top_TP[domidx_j])
					{
					m_score_top_TP[domidx_j] = Score;
					m_domidx_top_TP[domidx_j] = domidx_i;
					}
				}
			}
		else
			{
			++nf;
			if (m_score_top_FP[domidx_i] == FLT_MAX ||
				Score > m_score_top_FP[domidx_i])
				{
				m_score_top_FP[domidx_i] = Score;
				m_domidx_top_FP[domidx_i] = domidx_j;
				}
			if (triangle)
				{
				if (m_score_top_FP[domidx_j] == FLT_MAX ||
					Score > m_score_top_FP[domidx_j])
					{
					m_score_top_FP[domidx_j] = Score;
					m_domidx_top_FP[domidx_j] = domidx_i;
					}
				}
			}
		}
	float EPQ = 2*float(nf)/ndom;
	float Sens = 2*float(nt)/m_look->m_NT;

	if (m_SEPQ0_1 == FLT_MAX)
		{
		if (EPQ <= 0.1)	
			m_SEPQ0_1 = Sens;
		else
			m_SEPQ0_1 = 0;
		}

	if (m_SEPQ1 == FLT_MAX)
		{
		if (EPQ <= 1)	
			m_SEPQ1 = Sens;
		else
			m_SEPQ1 = 0;
		}

	if (m_SEPQ10 == FLT_MAX)
		{
		if (EPQ <= 0.1)	
			m_SEPQ10 = Sens;
		else
			m_SEPQ10 = 0;
		}

	m_CVESum3 = m_SEPQ0_1*2 + m_SEPQ1*3/2 + m_SEPQ10;
	//m_TopSum3 = m_top_SEPQ0_001*2 + m_top_SEPQ0_01*3/2 + m_top_SEPQ0_1;

	if (Msg != "noshow")
		{
		if (Msg != "")
			ProgressLog("%s ", Msg.c_str());
		ProgressLog("SEPQ0.1=%.3f", m_SEPQ0_1);
		ProgressLog(" SEPQ1=%.3f", m_SEPQ1);
		ProgressLog(" SEPQ10=%.3f", m_SEPQ10);
		ProgressLog(" Sum3=%.3f", m_CVESum3);
		ProgressLog(" %s", m_look->get_truthstr());
		ProgressLog("\n");

		//ProgressLog("TOPQ0.001=%.3f", m_top_SEPQ0_001);
		//ProgressLog(" TOPQ0.01=%.3f", m_top_SEPQ0_01);
		//ProgressLog(" TOPQ0.1=%.3f", m_top_SEPQ0_1);
		//ProgressLog(" Top3=%.3f", m_TopSum3);
		//ProgressLog(" %s", m_look->get_truthstr());
		//ProgressLog("\n");
		}

	asserta(!optset_top3);
	if (opt(top3))
		return m_TopSum3;
	else
		return m_CVESum3;
	}

void FastBench::ReadHits(
	const string &FN,
	uint qidx,
	uint tidx,
	uint scoreidx)
	{
	asserta(qidx > 0);
	asserta(tidx > 0);
	asserta(scoreidx > 0);
	--qidx;
	--tidx;
	--scoreidx;
	const uint maxidx = max(max(qidx, tidx), scoreidx);

#if SAVE_NOT_IN_DOPE
	FILE *fnid = CreateStdioFile("../tmp/tpnotindope.tmp");
#endif

	m_SeqCount = m_look->get_ndom();
	Alloc();
	const uint K = m_PairCount;
	for (uint i = 0; i < K; ++i)
		m_Scores[i] = FLT_MAX;
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
		uint qidx = m_look->get_domidx(q, true);
		uint tidx = m_look->get_domidx(t, true);
		if (qidx == UINT_MAX)
			missing.insert(q);
		if (tidx == UINT_MAX)
			missing.insert(t);
		if (qidx == UINT_MAX || tidx == UINT_MAX)
			continue;
		if (qidx == tidx)
			continue;
		uint k = triangle_ij_to_k2(qidx, tidx, m_SeqCount);
		assert(k < K);
		if (m_Scores[k] == FLT_MAX)
			{
			++n;
			bool tp = m_look->is_tp_k(k);
			if (tp) ++ntp;
			if (m_dope)
				{
				bool is_in_dope = in_dope(k);
				if (is_in_dope)
					{
					if (tp)
						++ntp_dope;
					else
						++nfp_dope;
					}
				else
					{
#if SAVE_NOT_IN_DOPE
					if (tp)
						{
						fputs(line.c_str(), fnid);
						fputc('\n', fnid);
						}
#endif
					continue;
					}
				}
			const float score = (float) StrToFloat(flds[scoreidx]);
			m_Scores[k] = score;
			}
		}
	Progress("Hits 100.00%%\n");
	ProgressLog("%u hits (%.3g%% of triangle), %u TPs\n",
		n, GetPct(n, K), ntp);
	if (m_dope)
		ProgressLog("%u TPs, %u FPs in dope\n", ntp_dope, nfp_dope);
	size_t nmiss = missing.size();
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
#if SAVE_NOT_IN_DOPE
	CloseStdioFile(fnid);
#endif
	}

void FastBench::WriteBits(const string &FN) const
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);
	uint K = triangle_get_K(m_SeqCount);
	WriteStdioFile(f, m_Scores, K*sizeof(m_Scores[0]));
	CloseStdioFile(f);
	}

void FastBench::ReadBits(const string &FN)
	{
	asserta(m_look);
	Alloc();
	FILE *f = OpenStdioFile(FN);
	uint K = triangle_get_K(m_SeqCount);
	ReadStdioFile(f, m_Scores, K*sizeof(m_Scores[0]));
	CloseStdioFile(f);
	}

void FastBench::WriteHits(const string &FN, bool IncludeSelf,
	bool UpperTriangleOnly) const
	{
	if (FN == "")
		return;
	asserta(m_ScoreOrder != 0);

	FILE *f = CreateStdioFile(FN);
	uint K = triangle_get_K(m_SeqCount);
	const vector<string> &labels = m_look->m_doms;
	for (uint k = 0; k < K; ++k)
		{
		ProgressStep(k, K, "Writing %s", FN.c_str());
		uint HitIdx = m_ScoreOrder[k];
		uint i, j;
		triangle_k_to_ij(HitIdx, m_SeqCount, i, j);
		if (i == j && !IncludeSelf)
			continue;

		float Score = m_Scores[HitIdx];
		if (Score == FLT_MAX) continue;

		fprintf(f, "%.3g", Score);
		fprintf(f, "\t%s", labels[i].c_str());
		fprintf(f, "\t%s", labels[j].c_str());
		fprintf(f, "\n");

		if (!UpperTriangleOnly)
			{
			fprintf(f, "%.3g", Score);
			fprintf(f, "\t%s", labels[j].c_str());
			fprintf(f, "\t%s", labels[i].c_str());
			fprintf(f, "\n");
			}
		}
	CloseStdioFile(f);
	}

void FastBench::ClearHitsAndResults()
	{
	m_CVESum3 = FLT_MAX;
	m_TopSum3 = FLT_MAX;
	m_SEPQ0_1 = FLT_MAX;
	m_SEPQ1 = FLT_MAX;
	m_SEPQ10 = FLT_MAX;
	m_top_SEPQ0_001 = FLT_MAX;
	m_top_SEPQ0_01 = FLT_MAX;
	m_top_SEPQ0_1 = FLT_MAX;
	const uint ndom = m_look->get_ndom();
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		m_score_top_TP[domidx] = FLT_MAX;
		m_score_top_FP[domidx] = FLT_MAX;
		m_domidx_top_TP[domidx] = UINT_MAX;
		m_domidx_top_FP[domidx] = UINT_MAX;
		}
	SubclassClearHitsAndResults();
	}

void FastBench::SetLookupFromLabels()
	{
	if (m_look == 0) m_look = new lookup;
	m_look->from_labels(m_Labels);
	}

void FastBench::ReadLookup(const string &argFN)
	{
	const string &FN = (argFN == "" ? "../data/scop40x.lookup" : argFN);
	m_scores_are_evalues = opt(scores_are_evalues);
	if (m_look == 0) m_look = new lookup;
	m_look->from_tsv(FN);
	m_SeqCount = m_look->get_ndom();
	m_PairCount = m_look->get_pair_count_upper_triangle_with_diagonal();
	m_Labels = m_look->m_doms;
	}

void FastBench::log_dope_ks() const
	{
	assert(m_look);
	const uint ndom = m_look->get_ndom();
	uint ntp = 0;
	uint nfp = 0;
	for (uint i = 0; i < m_dope_nhit; ++i)
		{
		ProgressStep(i, m_dope_nhit, "Logging dope hits");
		uint k = m_dope_ks[i];
		uint domidx_q, domidx_t;
		triangle_k_to_ij(k, ndom, domidx_q, domidx_t);
		string label_q, label_t;
		m_look->get_dom_scopid(domidx_q, label_q);
		m_look->get_dom_scopid(domidx_t, label_t);
		bool tp = m_look->is_tp_ij(domidx_q, domidx_t);
		if (tp) ntp++ ; else nfp++;
		Log("%s\t%s\t%s\n", label_q.c_str(), label_t.c_str(), tp ? "T" : "F");
		}
	ProgressLog("%u TP, %u FP\n", ntp, nfp);
	}

void FastBench::ReadDope(const string &FN)
	{
	uint8_t *read_bitdope(const string &fn, uint32_t &ndom, uint32_t &nhit);

	uint32_t ndom;
	m_dope = read_bitdope(FN, ndom, m_dope_nhit);
	asserta(ndom == m_look->get_ndom());
	m_dope_ks = myalloc(uint32_t, m_dope_nhit);
	const uint K = triangle_get_K(ndom);

	uint32_t bytes = (K + 7)/8;
	uint nhit = 0;
	for (uint i = 0; i < bytes; ++i)
		{
		uint8_t b = m_dope[i];
		for (uint j = 0; j < 8; ++j)
			{
			if (b & (1 << j))
				{
				uint k = i*8 + j;
				uint domidx_i, domidx_j;
				triangle_k_to_ij(k, ndom, domidx_i, domidx_j);
				assert(k < K);
				m_dope_ks[nhit++] = k;
				}
			}
		}
	asserta(nhit == m_dope_nhit);
	}

void cmd_fast_bench_hits()
	{
	asserta(!optset_output);
	const string &hitsfn = g_Arg1;
	const string lookupfn =
		(optset_lookup ? opt(lookup) : "../data/scop40x.lookup");

	FastBench FB;
	FB.ReadLookup(lookupfn);
	if (optset_dope)
		FB.ReadDope(opt(dope));
	uint qidx = 1;
	uint tidx = 2;
	uint scoreidx = 3;
	if (optset_qfield) qidx = opt(qfield);
	if (optset_tfield) tidx = opt(tfield);
	if (optset_scorefield) scoreidx = opt(scorefield);
	FB.ReadHits(hitsfn, qidx, tidx, scoreidx);
	FB.SetScoreOrder();
	FB.Bench();
	FB.WriteTopHits(opt(output2));
	}

void cmd_fast_bench_bits()
	{
	asserta(optset_lookup);
	asserta(!optset_output);
	const string &bitsfn = g_Arg1;

	FastBench FB;
	FB.ReadLookup(opt(lookup));
	FB.ReadBits(bitsfn);
	FB.SetScoreOrder();
	FB.Bench();
	}

void cmd_fb_hits2bits()
	{
	asserta(!optset_scorefirst);
	asserta(optset_lookup);
	asserta(optset_output);
	const string &hitsfn = g_Arg1;
	const string &outputfn = opt(output);

	uint qidx = 1;
	uint tidx = 2;
	uint scoreidx = 3;
	if (optset_qfield) qidx = opt(qfield);
	if (optset_tfield) tidx = opt(tfield);
	if (optset_scorefield) scoreidx = opt(scorefield);

	FastBench FB;
	FB.m_scores_are_evalues = opt(scores_are_evalues);
	FB.ReadLookup(opt(lookup));
	FB.ReadHits(hitsfn, qidx, tidx, scoreidx);
	FB.SetScoreOrder();
	FB.Bench();

	ProgressLog("Write %s\n", outputfn.c_str());
	FB.WriteBits(outputfn);
	}
