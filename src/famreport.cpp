#include "myutils.h"
#include "sort.h"
#include "scop40bench.h"

void SCOP40Bench::GetFamStats(uint FamIdx, 
  uint &HitIdxFirstFP, uint &TPsAtScoreFirstFP) const
	{
	HitIdxFirstFP = UINT_MAX;
	TPsAtScoreFirstFP = 0;
	const vector<uint> &DomIdxs = m_FamIdxToDomIdxs[FamIdx];
	float ScoreFirstFP = GetVeryBadScore();
	for (uint i = 0; i < SIZE(DomIdxs); ++i)
		{
		uint DomIdx = DomIdxs[i];
		float Score = m_DomIdxToScoreFirstFP[DomIdx];
		if (ScoreIsBetter(Score, ScoreFirstFP))
			{
			ScoreFirstFP = Score;
			HitIdxFirstFP = m_DomIdxToHitIdxFirstFP[DomIdx];
			}
		}

	for (uint i = 0; i < SIZE(DomIdxs); ++i)
		{
		uint DomIdx = DomIdxs[i];
		const vector<uint> &HitIdxs = m_DomIdxToHitIdxs[DomIdx];
		for (uint j = 0; j < SIZE(HitIdxs); ++j)
			{
			uint HitIdx = HitIdxs[j];
			uint DomIdx2 = m_DomIdx2s[HitIdx];
			if (m_DomIdxToFamIdx[DomIdx2] == FamIdx)
				{
				float Score = m_Scores[HitIdx];
				if (ScoreIsBetter(Score, ScoreFirstFP))
					++TPsAtScoreFirstFP;
				}
			}
		}
	}

void SCOP40Bench::FamReport(const string &FileName)
	{
	if (FileName.empty())
		return;
	FILE *f = CreateStdioFile(FileName);
	setbuf(f, 0);
	const uint FamCount = GetFamCount();
	SetFamSizes();
	vector<uint> Order(FamCount);
	QuickSortOrderDesc(m_FamSizes.data(), FamCount, Order.data());
	uint MinSize = 10;
	if (optset_minsize)
		MinSize = opt_minsize;
	double SumSens = 0;
	for (uint k = 0; k < FamCount; ++k)
		{
		ProgressStep(k, FamCount, "Reporting");
		uint FamIdx = Order[k];
		double MeanLength = GetMeanLength(FamIdx);
		uint Size = m_FamSizes[FamIdx];
		uint HitIdx_FirstFP, TPs_FirstFP;
		GetFamStats(FamIdx, HitIdx_FirstFP, TPs_FirstFP);
		if (Size < 10)
			break;
		uint DomIdx1_FirstFP = m_DomIdx1s[HitIdx_FirstFP];
		uint DomIdx2_FirstFP = m_DomIdx2s[HitIdx_FirstFP];
		double dSens = double(TPs_FirstFP)/Size;
		SumSens += dSens;
		fprintf(f, "%s", m_Fams[FamIdx].c_str());
		fprintf(f, "\t%.1f", MeanLength);
		fprintf(f, "\t%u", Size);
		fprintf(f, "\t%.3g", dSens);
		if (HitIdx_FirstFP == UINT_MAX)
			fprintf(f, "\t.\t.\t.");
		else
			{
			uint DomIdx1 = m_DomIdx1s[HitIdx_FirstFP];
			uint DomIdx2 = m_DomIdx2s[HitIdx_FirstFP];
			fprintf(f, "\t%s", m_Doms[DomIdx1].c_str());
			fprintf(f, "\t%s", m_Doms[DomIdx2].c_str());
			fprintf(f, "\t%.3g", m_Scores[HitIdx_FirstFP]);
			}
		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	ProgressLog("SumSens = %.1f\n", SumSens);
	}

void cmd_famreport()
	{
	//asserta(optset_labels);
	asserta(optset_input);
	SCOP40Bench SB;
	SB.ReadBit(g_Arg1);
	//SB.LoadLabels(opt_labels);
	vector<uint> SavedDomIdxToFamIdx = SB.m_DomIdxToFamIdx;
	SB.m_DomIdxToFamIdx.clear();
	SB.ReadChains(opt_input);
	asserta(SB.m_DomIdxToFamIdx == SavedDomIdxToFamIdx);
	SB.SetFamSizes();
	SB.ScanDomHits();
	SB.SetFamIdxToDomIdxs();
	SB.SetDomIdxToHitIdxs();
	SB.FamReport(opt_report);
	}
