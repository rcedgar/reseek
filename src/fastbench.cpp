#include "myutils.h"
#include "triangle.h"
#include "fastbench.h"
#include "sort.h"

void FastBench::Alloc()
	{
	m_SeqCount = SIZE(m_Labels);
	if (m_Scores != 0)
		myfree(m_Scores);
	m_PairCount = triangle_get_K(m_SeqCount) + 1;
	m_Scores = myalloc(float, m_PairCount);
	}

void FastBench::AppendLabel(const string &Label)
	{
	m_Labels.push_back(Label);
	}

void FastBench::AppendHit(uint i, uint j, float Score)
	{
	uint k = triangle_ij_to_k(i, j, m_SeqCount);
	m_Scores[k] = Score;
	SubclassAppendHit(i, j, Score);
	}

void FastBench::SetLookupFromLabels()
	{
	m_SFIdxToSize.clear();
	m_SFIdxToSize.resize(2000, 0);

	const uint N = SIZE(m_Labels);
	for (uint LabelIdx = 0; LabelIdx < N; ++LabelIdx)
		{
		const string &Label = m_Labels[LabelIdx];
		vector<string> Fields, Fields2;
		Split(Label, Fields, '/');
		asserta(SIZE(Fields) == 2);
		const string &Dom = Fields[0];
		const string &ScopId = Fields[1];
		Split(ScopId, Fields2, '.');
		asserta(SIZE(Fields2) == 3 || SIZE(Fields2) == 4);
		const string &SF = Fields2[0] + "." + Fields2[1] + "." + Fields2[2];
		AddDom(Dom, SF, LabelIdx);
		}

	m_NT = 0;
	m_NF = 0;
	const uint SFCount = SIZE(m_SFs);
	for (uint SFIdx = 0; SFIdx < SFCount; ++SFIdx)
		{
		uint Size = m_SFIdxToSize[SFIdx];
		asserta(Size > 0);
		m_NT += (Size*(Size - 1))/2;
		}
	m_NT *= 2;
	uint DomCount = SIZE(m_Labels);
	m_NF = DomCount*(DomCount-1) - m_NT;
	}


void FastBench::SetScoreOrder()
	{
	uint K = triangle_get_K(m_SeqCount);
	if (m_ScoreOrder == 0)
		myfree(m_ScoreOrder);
	m_ScoreOrder = myalloc(uint, K);
	QuickSortOrderDesc(m_Scores, K, m_ScoreOrder);
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
			float Sens = 2*float(nt)/m_NT;
			if (SEPQ0_1 == FLT_MAX && EPQ >= 0.1) SEPQ0_1 = Sens;
			if (SEPQ1 == FLT_MAX   && EPQ >= 1)   SEPQ1   = Sens;
			if (SEPQ10 == FLT_MAX  && EPQ >= 10)  SEPQ10  = Sens;
			LastScore = Score;
			}
		uint SFIdx_i = m_LabelIdxToSFIdx[LabelIdx_i];
		uint SFIdx_j = m_LabelIdxToSFIdx[LabelIdx_j];
		if (SFIdx_i == SFIdx_j)
			++nt;
		else
			++nf;
		}
	float EPQ = 2*float(nf)/m_SeqCount;
	float Sens = 2*float(nt)/m_NT;
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

void FastBench::AddDom(
	const string &Dom, const string &SF, uint LabelIdx)
	{
	uint SFIdx = UINT_MAX;
	if (m_SFToIdx.find(SF) == m_SFToIdx.end())
		{
		SFIdx = SIZE(m_SFs);
		m_SFs.push_back(SF);
		m_SFToIdx[SF] = SFIdx;
		}
	else
		SFIdx = m_SFToIdx[SF];
 
	m_LabelIdxToSFIdx.push_back(SFIdx);

	asserta(SFIdx < SIZE(m_SFIdxToSize));
	m_SFIdxToSize[SFIdx] += 1;
	}


void FastBench::WriteHits(const string &FN) const
	{
	if (FN == "")
		return;
	asserta(m_ScoreOrder != 0);

	FILE *f = CreateStdioFile(FN);
	uint K = triangle_get_K(m_SeqCount);
	for (uint k = 0; k < K; ++k)
		{
		ProgressStep(k, K, "Writing %s", FN.c_str());
		uint HitIdx = m_ScoreOrder[k];
		uint i, j;
		triangle_k_to_ij(HitIdx, m_SeqCount, i, j);
		if (i == j)
			continue;

		fprintf(f, "%.3g", m_Scores[HitIdx]);
		fprintf(f, "\t%s", m_Labels[i].c_str());
		fprintf(f, "\t%s", m_Labels[j].c_str());
		fprintf(f, "\n");

		fprintf(f, "%.3g", m_Scores[HitIdx]);
		fprintf(f, "\t%s", m_Labels[j].c_str());
		fprintf(f, "\t%s", m_Labels[i].c_str());
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
