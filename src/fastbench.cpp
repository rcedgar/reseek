#include "myutils.h"
#include "triangle.h"
#include "fastbench.h"
#include "sort.h"
#include <algorithm>
#include <unordered_set>

#define SAVE_NOT_IN_DOPE	0

#if PARALLEL_SORT

#include "parallel_sort.h"

void FastBench::SetScoreOrder_Parallel()
{
	asserta(m_Scores);
	const uint K = triangle_get_K(m_ndom);

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

#if PARALLEL_SORT
	if (m_Scores == 0)
		{
		m_Scores = myalloc(float, npair);
		for (uint i = 0; i < npair; ++i) m_Scores[i] = FLT_MAX;
		m_npair = npair;
		}
	else
		asserta(m_npair == npair);
#else
	m_PairCount = npair;
	myfree(m_Scores);
	myfree(m_ScoreOrder);
	m_Scores = myalloc(float, npair);
#endif
	}

void FastBench::AppendHit(uint i, uint j, float Score)
	{
	uint k = triangle_ij_to_k(i, j, m_ndom);
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
	uint K = triangle_get_K(m_ndom);
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

bool FastBench::IsTprHoldout(uint domidx_q, uint domidx_t) const
	{
	assert(m_look);
	if (m_look->m_LT == LT_TOP_SF)
		return m_look->same_fam_ij(domidx_q, domidx_t);
	if (m_look->m_LT == LT_TOP_FOLD)
		return m_look->same_sf_ij(domidx_q, domidx_t);
	return false;
	}

bool FastBench::IsFprHoldout(uint domidx_q, uint domidx_t) const
	{
	assert(m_look);
	if (m_look->m_LT == LT_TOP_SF)
		return m_look->same_sf_ij(domidx_q, domidx_t);
	if (m_look->m_LT == LT_TOP_FOLD)
		return m_look->same_fold_ij(domidx_q, domidx_t);
	return false;
	}

bool FastBench::BetterScore(float score1, float score2) const
	{
	if (score2 == FLT_MAX)
		return score1 != FLT_MAX;
	if (score1 == FLT_MAX)
		return false;
	if (m_scores_are_evalues)
		return score1 < score2;
	return score1 > score2;
	}

bool FastBench::PassesThreshold(float score, float tau) const
	{
	if (score == FLT_MAX)
		return false;
	if (m_scores_are_evalues)
		return score <= tau;
	return score >= tau;
	}

double FastBench::BenchTop(const string &Msg)
	{
	asserta(m_look);
	asserta(m_Scores);
	const uint ndom = m_look->get_ndom();
	const uint K = triangle_get_K(ndom);

	const bool top_sf = (m_look->m_LT == LT_TOP_SF);
	const bool top_fold = (m_look->m_LT == LT_TOP_FOLD);
	asserta(top_sf || top_fold);

	vector<float> TprBestTp(ndom, FLT_MAX);
	vector<float> TprBestFp(ndom, FLT_MAX);
	vector<float> FprBestAny(ndom, FLT_MAX);

	for (uint k = 0; k < K; ++k)
		{
		float Score = m_Scores[k];
		if (Score == FLT_MAX)
			continue;
		uint domidx_i, domidx_j;
		triangle_k_to_ij(k, ndom, domidx_i, domidx_j);
		if (domidx_i == domidx_j)
			continue;

		auto UpdateQuery = [&](uint domidx_q, uint domidx_t)
			{
			if (!IsTprHoldout(domidx_q, domidx_t))
				{
				if (IsTP(domidx_q, domidx_t))
					{
					if (BetterScore(Score, TprBestTp[domidx_q]))
						TprBestTp[domidx_q] = Score;
					}
				else if (BetterScore(Score, TprBestFp[domidx_q]))
					{
					TprBestFp[domidx_q] = Score;
					}
				}
			if (!IsFprHoldout(domidx_q, domidx_t))
				{
				if (BetterScore(Score, FprBestAny[domidx_q]))
					FprBestAny[domidx_q] = Score;
				}
			};

		UpdateQuery(domidx_i, domidx_j);
		UpdateQuery(domidx_j, domidx_i);
		}

	vector<bool> InQ(ndom, false);
	if (top_sf)
		{
		const uint nsf = uint(m_look->m_sfs.size());
		vector<unordered_set<uint> > sf_fams(nsf);
		for (uint domidx = 0; domidx < ndom; ++domidx)
			{
			uint sfidx = m_look->m_domidx2sfidx[domidx];
			uint famidx = m_look->m_domidx2famidx[domidx];
			sf_fams[sfidx].insert(famidx);
			}
		for (uint domidx = 0; domidx < ndom; ++domidx)
			{
			uint sfidx = m_look->m_domidx2sfidx[domidx];
			if (sf_fams[sfidx].size() >= 2)
				InQ[domidx] = true;
			}
		}
	else
		{
		const uint nfold = uint(m_look->m_folds.size());
		vector<unordered_set<uint> > fold_sfs(nfold);
		for (uint domidx = 0; domidx < ndom; ++domidx)
			{
			uint foldidx = m_look->m_domidx2foldidx[domidx];
			uint sfidx = m_look->m_domidx2sfidx[domidx];
			fold_sfs[foldidx].insert(sfidx);
			}
		for (uint domidx = 0; domidx < ndom; ++domidx)
			{
			uint foldidx = m_look->m_domidx2foldidx[domidx];
			if (fold_sfs[foldidx].size() >= 2)
				InQ[domidx] = true;
			}
		}

	uint n_q = 0;
	for (uint domidx = 0; domidx < ndom; ++domidx)
		if (InQ[domidx])
			++n_q;
	asserta(n_q > 0);

	auto GetBestTprAny = [&](uint domidx) -> float
		{
		const float best_tp = TprBestTp[domidx];
		const float best_fp = TprBestFp[domidx];
		if (best_tp == FLT_MAX && best_fp == FLT_MAX)
			return FLT_MAX;
		if (best_tp == FLT_MAX)
			return best_fp;
		if (best_fp == FLT_MAX)
			return best_tp;
		if (BetterScore(best_tp, best_fp))
			return best_tp;
		return best_fp;
		};

	auto TopIsTp = [&](uint domidx) -> bool
		{
		const float best_tp = TprBestTp[domidx];
		const float best_fp = TprBestFp[domidx];
		if (best_tp == FLT_MAX)
			return false;
		if (best_fp == FLT_MAX)
			return true;
		return BetterScore(best_tp, best_fp) || best_tp == best_fp;
		};

	vector<float> taus;
	taus.reserve(n_q*2);
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		if (!InQ[domidx])
			continue;
		const float best_tpr = GetBestTprAny(domidx);
		if (best_tpr != FLT_MAX)
			taus.push_back(best_tpr);
		if (FprBestAny[domidx] != FLT_MAX)
			taus.push_back(FprBestAny[domidx]);
		}
	asserta(taus.size() > 0);
	QuickSortInPlace(taus.data(), uint(taus.size()));
	taus.erase(unique(taus.begin(), taus.end()), taus.end());
	if (!m_scores_are_evalues)
		reverse(taus.begin(), taus.end());

	m_Sum3 = FLT_MAX;
	m_SEPQtopA = FLT_MAX;
	m_SEPQtopB = FLT_MAX;
	m_SEPQtopC = FLT_MAX;

	const float FprThreshA = 0.001f;
	const float FprThreshB = 0.01f;
	const float FprThreshC = 0.1f;

	FILE *fRoc = 0;
	if (optset_roc)
		fRoc = CreateStdioFile(opt(roc));

	auto WriteCvePoint = [&](float Tpr, float Fpr, float Score)
		{
		if (fRoc == 0)
			return;
		fprintf(fRoc, "%.6g", Fpr);
		fprintf(fRoc, "\t%.6g", Tpr);
		fprintf(fRoc, "\t%.6g", Score);
		fprintf(fRoc, "\n");
		};

	auto UpdateSepq = [&](float Fpr, float Tpr)
		{
		if (Fpr >= FprThreshA && m_SEPQtopA == FLT_MAX)
			m_SEPQtopA = Tpr;
		if (Fpr >= FprThreshB && m_SEPQtopB == FLT_MAX)
			m_SEPQtopB = Tpr;
		if (Fpr >= FprThreshC && m_SEPQtopC == FLT_MAX)
			m_SEPQtopC = Tpr;
		};

	auto ClassifyAtTau = [&](float tau, uint &n_tp, uint &n_fn,
		uint &n_fp, uint &n_tn)
		{
		n_tp = n_fn = n_fp = n_tn = 0;
		for (uint domidx = 0; domidx < ndom; ++domidx)
			{
			if (!InQ[domidx])
				continue;
			const float best_tpr = GetBestTprAny(domidx);
			if (!PassesThreshold(best_tpr, tau))
				++n_fn;
			else if (TopIsTp(domidx))
				++n_tp;

			if (!PassesThreshold(FprBestAny[domidx], tau))
				++n_tn;
			else
				++n_fp;
			}
		};

	if (fRoc)
		fprintf(fRoc, "FPR\tTPR\t%s\n",
			m_scores_are_evalues ? "evalue" : "score");
	//WriteCvePoint(0.0f, 0.0f,
	//	m_scores_are_evalues ? 0.0f : FLT_MAX);

	float final_tpr = 0.0f;
	float final_fpr = 0.0f;
	for (uint ti = 0; ti < uint(taus.size()); ++ti)
		{
		const float tau = taus[ti];
		uint n_tp, n_fn, n_fp, n_tn;
		ClassifyAtTau(tau, n_tp, n_fn, n_fp, n_tn);
		const float tpr = float(n_tp)/float(n_tp + n_fn);
		const float fpr = float(n_fp)/float(n_fp + n_tn);
		WriteCvePoint(tpr, fpr, tau);
		UpdateSepq(fpr, tpr);
		final_tpr = tpr;
		final_fpr = fpr;
		}

	{
	const float loose_tau = m_scores_are_evalues ? FLT_MAX : -FLT_MAX;
	uint n_tp, n_fn, n_fp, n_tn;
	ClassifyAtTau(loose_tau, n_tp, n_fn, n_fp, n_tn);
	final_tpr = float(n_tp)/float(n_tp + n_fn);
	final_fpr = float(n_fp)/float(n_fp + n_tn);
	WriteCvePoint(final_tpr, final_fpr, loose_tau);
	UpdateSepq(final_fpr, final_tpr);
	}

	if (fRoc)
		CloseStdioFile(fRoc);

	if (m_SEPQtopA == FLT_MAX) m_SEPQtopA = final_tpr;
	if (m_SEPQtopB == FLT_MAX) m_SEPQtopB = final_tpr;
	if (m_SEPQtopC == FLT_MAX) m_SEPQtopC = final_tpr;

	m_Sum3 = m_SEPQtopA*2.0f + m_SEPQtopB*1.5f + m_SEPQtopC;

	if (Msg != "noshow")
		{
		if (Msg != "")
			ProgressLog("%s ", Msg.c_str());
		ProgressLog("TEPQ%.3g=%.3f", FprThreshA, m_SEPQtopA);
		ProgressLog(" TEPQ%.3g=%.3f", FprThreshB, m_SEPQtopB);
		ProgressLog(" TEPQ%.3g=%.3f", FprThreshC, m_SEPQtopC);
		ProgressLog(" Top3=%.3f", m_Sum3);
		ProgressLog(" %s", m_look->get_truthstr());
		if (m_name != "")
			ProgressLog(" %s", m_name.c_str());
		ProgressLog("\n");
		}

	asserta(!optset_top3);
	return m_Sum3;
	}

double FastBench::Bench(const string &Msg)
	{
	asserta(m_look);
	if (m_look->m_LT == LT_TOP_SF || m_look->m_LT == LT_TOP_FOLD)
		return BenchTop(Msg);
	return BenchPairCVE(Msg);
	}

// This is not ROC, it is Coverage vs. Error (CVE)
double FastBench::BenchPairCVE(const string &Msg)
	{
	asserta(m_ScoreOrder != 0);
	const uint ndom = m_look->get_ndom();
	uint K = triangle_get_K(ndom);
	uint nt = 0;
	uint nf = 0;

	m_Sum3 = FLT_MAX;
	m_SEPQ0_1 = FLT_MAX;
	m_SEPQ1 = FLT_MAX;
	m_SEPQ10 = FLT_MAX;

	float LastScore = m_scores_are_evalues ? -9e9f : FLT_MAX;

	m_SEPQ0_1 = FLT_MAX;
	m_SEPQ1 = FLT_MAX;
	m_SEPQ10 = FLT_MAX;

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

		// must check when score changes
		// to avoid sort artifacts!
		if (Score != LastScore)
			{
			if (m_scores_are_evalues)
				asserta(Score > LastScore);
			else
				asserta(Score < LastScore);
			float EPQ = 2*float(nf)/ndom;
			float Sens = 2*float(nt)/m_look->m_NT;
			if (EPQ <= 0.1) m_SEPQ0_1 = Sens;
			if (EPQ <= 1)   m_SEPQ1   = Sens;
			if (EPQ <= 10)  m_SEPQ10  = Sens;

			LastScore = Score;
			}
		if (IsIgnored(domidx_i, domidx_j))
			continue;

		if (IsTP(domidx_i, domidx_j))
			++nt;
		else
			++nf;
		}
	float EPQ = 2*float(nf)/ndom;
	float Sens = 2*float(nt)/m_look->m_NT;

	m_Sum3 = m_SEPQ0_1*2 + m_SEPQ1*3/2 + m_SEPQ10;

	if (Msg != "noshow")
		{
		if (Msg != "")
			ProgressLog("%s ", Msg.c_str());
		ProgressLog("SEPQ0.1=%.3f", m_SEPQ0_1);
		ProgressLog(" SEPQ1=%.3f", m_SEPQ1);
		ProgressLog(" SEPQ10=%.3f", m_SEPQ10);
		ProgressLog(" Sum3=%.3f", m_Sum3);
		ProgressLog(" %s", m_look->get_truthstr());
		if (m_name != "")
			ProgressLog(" %s", m_name.c_str());
		ProgressLog("\n");
		}

	asserta(!optset_top3);
	return m_Sum3;
	}

double FastBench::calc_sffp()
	{
	asserta(m_look);
	asserta(m_Scores);
	const uint ndom = m_look->get_ndom();
	const uint K = triangle_get_K(ndom);
	const uint NT = m_look->m_NT;
	if (NT == 0)
		{
		m_SFFP = 0;
		return 0;
		}

	// Pass 1: find best FP score per query
	vector<float> best_fp(ndom, FLT_MAX);
	for (uint k = 0; k < K; ++k)
		{
		float Score = m_Scores[k];
		if (Score == FLT_MAX)
			continue;
		uint domidx_i, domidx_j;
		triangle_k_to_ij(k, ndom, domidx_i, domidx_j);
		if (domidx_i == domidx_j)
			continue;
		if (IsIgnored(domidx_i, domidx_j))
			continue;
		if (!IsTP(domidx_i, domidx_j))
			{
			if (BetterScore(Score, best_fp[domidx_i]))
				best_fp[domidx_i] = Score;
			if (BetterScore(Score, best_fp[domidx_j]))
				best_fp[domidx_j] = Score;
			}
		}

	// Pass 2: count TPs per query with score strictly better than best FP
	uint total_tp_above = 0;
	for (uint k = 0; k < K; ++k)
		{
		float Score = m_Scores[k];
		if (Score == FLT_MAX)
			continue;
		uint domidx_i, domidx_j;
		triangle_k_to_ij(k, ndom, domidx_i, domidx_j);
		if (domidx_i == domidx_j)
			continue;
		if (IsIgnored(domidx_i, domidx_j))
			continue;
		if (IsTP(domidx_i, domidx_j))
			{
			if (best_fp[domidx_i] == FLT_MAX || BetterScore(Score, best_fp[domidx_i]))
				++total_tp_above;
			if (best_fp[domidx_j] == FLT_MAX || BetterScore(Score, best_fp[domidx_j]))
				++total_tp_above;
			}
		}

	m_SFFP = float(total_tp_above) / float(NT);
	return m_SFFP;
	}

void FastBench::ReadHits(
	const string &FN,
	uint qidx,
	uint tidx,
	uint scoreidx)
	{
	GetStemName(FN, m_name);
	const uint maxidx = max(max(qidx, tidx), scoreidx);

#if SAVE_NOT_IN_DOPE
	FILE *fnid = CreateStdioFile("../tmp/tpnotindope.tmp");
#endif

	m_ndom = m_look->get_ndom();
	Alloc();
	const uint K = m_npair;
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
	uint linenr = 0;
	while (ReadLineStdioFile(f, line))
		{
		++linenr;
		if (line[0] == '#' || StartsWith(line, "q"))
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
		if (flds.size() < maxidx)
			Die("line %u not enough fields (%u) '%s'",
				linenr, uint(flds.size()), line.c_str());
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
		uint k = triangle_ij_to_k2(qidx, tidx, m_ndom);
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
	Progress("Hits 100.00%% (%s)\n", MemBytesToStr(n));
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
	uint K = triangle_get_K(m_ndom);
	WriteStdioFile(f, m_Scores, K*sizeof(m_Scores[0]));
	CloseStdioFile(f);
	}

void FastBench::ReadBits(const string &FN)
	{
	asserta(m_look);
	Alloc();
	FILE *f = OpenStdioFile(FN);
	uint K = triangle_get_K(m_ndom);
	ReadStdioFile(f, m_Scores, K*sizeof(m_Scores[0]));
	CloseStdioFile(f);
	}

void FastBench::WriteHits(const string &FN, bool IncludeSelf,
	bool UpperTriangleOnly, bool IncludeFam) const
	{
	if (FN == "")
		return;
	asserta(m_ScoreOrder != 0);

	FILE *f = CreateStdioFile(FN);
	uint K = triangle_get_K(m_ndom);
	const vector<string> &labels = m_look->m_doms;
	for (uint k = 0; k < K; ++k)
		{
		ProgressStep(k, K, "Writing %s include_self=%c triangle=%c fam=%c",
			FN.c_str(), tof(IncludeSelf), tof(UpperTriangleOnly), tof(IncludeFam));
		uint HitIdx = m_ScoreOrder[k];
		uint i, j;
		triangle_k_to_ij(HitIdx, m_ndom, i, j);
		if (i == j && !IncludeSelf)
			continue;

		float Score = m_Scores[HitIdx];
		if (Score == FLT_MAX) continue;

		fprintf(f, "%.3g", Score);
		fprintf(f, "\t%s", labels[i].c_str());
		if (IncludeFam)
			fprintf(f, "/%s", m_look->get_fam(i));
		fprintf(f, "\t%s", labels[j].c_str());
		if (IncludeFam)
			fprintf(f, "/%s", m_look->get_fam(j));
		fprintf(f, "\n");

		if (!UpperTriangleOnly)
			{
			fprintf(f, "%.3g", Score);
			fprintf(f, "\t%s", labels[j].c_str());
			if (IncludeFam)
				fprintf(f, "/%s", m_look->get_fam(j));
			fprintf(f, "\t%s", labels[i].c_str());
			if (IncludeFam)
				fprintf(f, "/%s", m_look->get_fam(i));
			fprintf(f, "\n");
			}
		}
	CloseStdioFile(f);
	}

void FastBench::ClearHitsAndResults()
	{
	m_Sum3 = FLT_MAX;
	m_SFFP = FLT_MAX;
	m_SEPQ0_1 = FLT_MAX;
	m_SEPQ1 = FLT_MAX;
	m_SEPQ10 = FLT_MAX;
	m_SEPQtopA = FLT_MAX;
	m_SEPQtopB = FLT_MAX;
	m_SEPQtopC = FLT_MAX;
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
	if (m_look == 0) m_look = new lookup;
	m_look->from_tsv(FN);
	m_ndom = m_look->get_ndom();
	m_npair = m_look->get_pair_count_upper_triangle_with_diagonal();
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

void guess_fields(
	const string &hitsfn,
	uint &qfi, uint &tfi, uint &sfi,
	bool &scores_are_evalues)
	{
	qfi = tfi = sfi = UINT_MAX;
	scores_are_evalues = false;

	FILE *f = OpenStdioFile(hitsfn);
	string line;
	vector<string> flds;
	uint n0 = 0;
	uint n1 = 0;
	uint n2 = 0;
	uint nlt1 = 0;
	for (uint i = 0; i < 100; ++i)
		{
		bool ok = ReadLineStdioFile(f, line);
		asserta(ok);
		if (line[0] == '#') continue;
		Split(line, flds, '\t');
		asserta(flds.size() >= 3);
		if (IsValidFloatStr(flds[0]))
			{
			++n0;
			if (StrToFloat(flds[0]) < 1)
				++nlt1;
			}
		if (IsValidFloatStr(flds[1]))
			{
			++n1;
			if (StrToFloat(flds[1]) < 1)
				++nlt1;
			}
		if (IsValidFloatStr(flds[2]))
			{
			++n2;
			if (StrToFloat(flds[2]) < 1)
				++nlt1;
			}
		}
	CloseStdioFile(f);
	if (n0 == 0 && n1 == 0 && n2 == 100)
		{
		qfi = 0;
		tfi = 1;
		sfi = 2;
		}
	else if (n0 == 100 && n1 == 0 && n2 == 0)
		{
		qfi = 1;
		tfi = 2;
		sfi = 0;
		}
	if (opt(scores_are_evalues))
		scores_are_evalues = true;
	else if (opt(scores_are_not_evalues))
		scores_are_evalues = false;
	else
		scores_are_evalues = (nlt1 > 10);
	if (qfi == UINT_MAX)
		Die("guess_fields");
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
			hitsfn, qidx, tidx, scoreidx, FB.m_scores_are_evalues);

	ProgressLog("%s(%u,%u,%u)\n",
		FB.m_scores_are_evalues ? "E-values" : "scores",
		qidx+1, tidx+1, scoreidx+1);

	FB.ReadHits(hitsfn, qidx, tidx, scoreidx);
	if (FB.m_look->m_LT != LT_TOP_SF && FB.m_look->m_LT != LT_TOP_FOLD)
		FB.SetScoreOrder();
	FB.Bench();
	if (opt(sffp))
		{
		double sffp = FB.calc_sffp();
		ProgressLog("SFFP=%.4f\n", sffp);
		}
	}

void cmd_fast_bench_bits()
	{
	asserta(optset_lookup);
	asserta(!optset_output);
	const string &bitsfn = g_Arg1;

	FastBench FB;
	FB.ReadLookup(opt(lookup));
	FB.ReadBits(bitsfn);
	if (FB.m_look->m_LT != LT_TOP_SF && FB.m_look->m_LT != LT_TOP_FOLD)
		FB.SetScoreOrder();
	FB.Bench();
	}
