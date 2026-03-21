#include "myutils.h"
#include "triangle.h"
#include "fastbench.h"
#include "sort.h"

void FastBench::Alloc()
	{
	asserta(m_look);
	const uint ndom = m_look->get_ndom();
	const uint npair = m_look->get_pair_count_upper_triangle_with_diagonal();
	if (m_Scores != 0)
		myfree(m_Scores);
	m_PairCount = triangle_get_K(ndom) + 1;
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
	uint K = triangle_get_K(m_SeqCount);
	if (m_ScoreOrder == 0)
		myfree(m_ScoreOrder);
	m_ScoreOrder = myalloc(uint, K);
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
	float LastScore = FLT_MAX;
	float SEPQ0_1 = FLT_MAX;
	float SEPQ1 = FLT_MAX;
	float SEPQ10 = FLT_MAX;
	for (uint k = 0; k < K; ++k)
		{
		uint HitIdx = m_ScoreOrder[k];
		uint LabelIdx_i, LabelIdx_j;
		triangle_k_to_ij(HitIdx, m_SeqCount, LabelIdx_i, LabelIdx_j);
		if (LabelIdx_i == LabelIdx_j)
			continue;
		float Score = m_Scores[HitIdx];
		if (Score != LastScore)
			{
			asserta(Score < LastScore);
			float EPQ = 2*float(nf)/m_SeqCount;
			float Sens = 2*float(nt)/m_look->m_NT;
			if (SEPQ0_1 == FLT_MAX && EPQ >= 0.1) SEPQ0_1 = Sens;
			if (SEPQ1 == FLT_MAX   && EPQ >= 1)   SEPQ1   = Sens;
			if (SEPQ10 == FLT_MAX  && EPQ >= 10)  SEPQ10  = Sens;
			LastScore = Score;
			}
		if (IsTP(LabelIdx_i, LabelIdx_j))
			++nt;
		else
			++nf;
		}
	float EPQ = 2*float(nf)/m_SeqCount;
	float Sens = 2*float(nt)/m_look->m_NT;
	if (SEPQ0_1 == FLT_MAX && EPQ >= 0.1) SEPQ0_1 = Sens;
	if (SEPQ1 == FLT_MAX   && EPQ >= 1)   SEPQ1   = Sens;
	if (SEPQ10 == FLT_MAX  && EPQ >= 10)  SEPQ10  = Sens;
	m_Sum3 = SEPQ0_1*2 + SEPQ1*3/2 + SEPQ10;

	if (Msg != "")
		ProgressLog("%s ", Msg.c_str());
	ProgressLog("SEPQ0.1=%.3f", SEPQ0_1);
	ProgressLog(" SEPQ1=%.3f", SEPQ1);
	ProgressLog(" SEPQ10=%.3f", SEPQ10);
	ProgressLog(" Sum3=%.3f", m_Sum3);
	ProgressLog("\n");
	}

void FastBench::WriteHits(const string &FN, bool IncludeSelf) const
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

		fprintf(f, "%.3g", m_Scores[HitIdx]);
		fprintf(f, "\t%s", labels[i].c_str());
		fprintf(f, "\t%s", labels[j].c_str());
		fprintf(f, "\n");

		fprintf(f, "%.3g", m_Scores[HitIdx]);
		fprintf(f, "\t%s", labels[j].c_str());
		fprintf(f, "\t%s", labels[i].c_str());
		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	}

void FastBench::ClearHitsAndResults()
	{
	myfree(m_Scores);
	m_Scores = 0;
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
	if (m_look == 0) m_look = new lookup;
	m_look->from_tsv(FN);
	}
