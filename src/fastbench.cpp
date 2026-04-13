#include "myutils.h"
#include "triangle.h"
#include "fastbench.h"
#include "sort.h"

#define SAVE_NOT_IN_DOPE	0

#if PARALLEL_SORT

#include "parallel_sort.h"

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
	myfree(m_Scores);
	myfree(m_ScoreOrder);
	m_PairCount = npair;
	m_Scores = myalloc(float, npair);
	}

void FastBench::AppendHit(uint i, uint j, float Score)
	{
	uint k = triangle_ij_to_k(i, j, m_SeqCount);
	m_Scores[k] = Score;
	SubclassAppendHit(i, j, Score);
	}

void FastBench::SetScoreOrder()
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

bool FastBench::IsTP(uint LabelIdx_i, uint LabelIdx_j) const
	{
	assert(m_look);
	return m_look->is_tp_ij(LabelIdx_i, LabelIdx_j);
	}

void FastBench::Bench(const string &Msg)
	{
	asserta(m_ScoreOrder != 0);
	uint K = triangle_get_K(m_SeqCount);
	uint nt = 0;
	uint nf = 0;
	float LastScore = m_scores_are_evalues ? -1 : FLT_MAX;
	m_SEPQ0_1 = FLT_MAX;
	m_SEPQ1 = FLT_MAX;
	m_SEPQ10 = FLT_MAX;
	for (uint k = 0; k < K; ++k)
		{
		uint HitIdx = m_ScoreOrder[k];
		uint LabelIdx_i, LabelIdx_j;
		triangle_k_to_ij(HitIdx, m_SeqCount, LabelIdx_i, LabelIdx_j);
		if (LabelIdx_i == LabelIdx_j)
			continue;
		float Score = m_Scores[HitIdx];
		if (Score == FLT_MAX) continue;
		if (Score != LastScore)
			{
			if (m_scores_are_evalues)
				asserta(Score > LastScore);
			else
				asserta(Score < LastScore);
			float EPQ = 2*float(nf)/m_SeqCount;
			float Sens = 2*float(nt)/m_look->m_NT;
			if (m_SEPQ0_1 == FLT_MAX && EPQ >= 0.1) m_SEPQ0_1 = Sens;
			if (m_SEPQ1 == FLT_MAX   && EPQ >= 1)   m_SEPQ1   = Sens;
			if (m_SEPQ10 == FLT_MAX  && EPQ >= 10)  m_SEPQ10  = Sens;
			LastScore = Score;
			}
		if (IsTP(LabelIdx_i, LabelIdx_j))
			++nt;
		else
			++nf;
		}
	float EPQ = 2*float(nf)/m_SeqCount;
	float Sens = 2*float(nt)/m_look->m_NT;
	if (m_SEPQ0_1 == FLT_MAX) m_SEPQ0_1 = Sens;
	if (m_SEPQ1 == FLT_MAX)   m_SEPQ1   = Sens;
	if (m_SEPQ10 == FLT_MAX)  m_SEPQ10  = Sens;
	m_Sum3 = m_SEPQ0_1*2 + m_SEPQ1*3/2 + m_SEPQ10;

	if (Msg != "noshow")
		{
		if (Msg != "")
			ProgressLog("%s ", Msg.c_str());
		ProgressLog("SEPQ0.1=%.3f", m_SEPQ0_1);
		ProgressLog(" SEPQ1=%.3f", m_SEPQ1);
		ProgressLog(" SEPQ10=%.3f", m_SEPQ10);
		ProgressLog(" Sum3=%.3f", m_Sum3);
		ProgressLog("\n");
		}
	}

void FastBench::ReadHits(
	const string &FN,
	uint qidx,
	uint tidx,
	uint scoreidx)
	{
#if SAVE_NOT_IN_DOPE
	FILE *fnid = CreateStdioFile("../tmp/tpnotindope.tmp");
#endif
	const uint maxidx = max(max(qidx, tidx), scoreidx);
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
	while (ReadLineStdioFile(f, line))
		{
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
		const string &t = flds[tidx];
		uint qidx = m_look->get_domidx(q);
		uint tidx = m_look->get_domidx(t);
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
	const vector<string> labels = m_look->m_doms;
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
	//myfree(m_Scores);
	//m_Scores = 0;
	m_Sum3 = FLT_MAX;
	SubclassClearHitsAndResults();
	}

void FastBench::SetLookupFromLabels()
	{
	if (m_look == 0) m_look = new lookup;
	m_look->from_labels(m_Labels);
	}

void FastBench::ReadLookup(const string &FN)
	{
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
	asserta(optset_lookup);
	asserta(!optset_output);
	const string &hitsfn = g_Arg1;

	FastBench FB;
	FB.ReadLookup(opt(lookup));
	if (optset_dope)
		FB.ReadDope(opt(dope));
	if (opt(scorefirst))
		FB.ReadHits(hitsfn, 1, 2, 0);
	else
		FB.ReadHits(hitsfn, 0, 1, 2);
	FB.SetScoreOrder();
	FB.Bench();
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

	uint qidx = 0;
	uint tidx = 1;
	uint scoreidx = 2;
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
