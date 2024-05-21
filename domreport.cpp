#include "myutils.h"
#include "sort.h"
#include "scop40bench.h"

void SCOP40Bench::DomReport(const string &FileName)
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
	uint SumSens = 0;
	for (uint k = 0; k < FamCount; ++k)
		{
		ProgressStep(k, FamCount, "Reporting");
		uint FamIdx = Order[k];
		vector<uint> &DomIdxs = m_FamIdxToDomIdxs[FamIdx];
		uint FamSize = m_FamSizes[FamIdx];
		asserta(SIZE(DomIdxs) == FamSize);
		for (uint i = 0; i < FamSize; ++i)
			{
			uint DomIdx = DomIdxs[i];
			const string &Label = m_Doms[DomIdx];
			float Score1FP = m_DomIdxToScoreFirstFP[DomIdx];
			uint Sens1FP = m_DomIdxToSens1FP[DomIdx];
			SumSens += Sens1FP;
			if (FamSize >= MinSize)
				{
				fprintf(f, "%s\t%u\t%.3g\n",
				  Label.c_str(), Sens1FP, Score1FP);
				}
			}
		}
	CloseStdioFile(f);
	ProgressLog("Dom SumSens1FP %u\n", SumSens);
	}

void cmd_domreport()
	{
	asserta(optset_input);
	SCOP40Bench SB;
	SB.ReadBit(g_Arg1);
	SB.m_DomIdxToFamIdx.clear();
	SB.ReadChains(opt_input);
	SB.SetFamSizes();
	SB.ScanDomHits();
	SB.SetFamIdxToDomIdxs();
	SB.SetDomIdxToHitIdxs();
	SB.DomReport(opt_report);
	}
