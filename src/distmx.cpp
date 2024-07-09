#include "myutils.h"
#include "dbsearcher.h"
#include "chainreader2.h"
#include "mx.h"
#include "pdbchain.h"
#include "xdpmem.h"
#include "cigar.h"
#include "timing.h"
#include "sort.h"
#include "distmxsearcher.h"

void DistMxSearcher::OnSetup()
	{
	asserta(m_fDistMx != 0);
	m_ChainCount = GetDBChainCount();
	fprintf(m_fDistMx, "distmx\t%u\n", m_ChainCount);
	for (uint ChainIndex = 0; ChainIndex < m_ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, m_ChainCount, "Write labels");
		const PDBChain &Chain = *m_DBChains[ChainIndex];
		const string &Label = Chain.m_Label;
		fprintf(m_fDistMx, "%u\t%s\n", ChainIndex, Label.c_str());
		}
	}

void DistMxSearcher::OnAln(DSSAligner &DA, bool Up)
	{
	if (!Up)
		return;
	asserta(m_fDistMx!= 0);
	float MeanTS = (DA.m_TestStatisticA + DA.m_TestStatisticB)/2;
	if (MeanTS < m_MinTS)
		return;
	m_MaxTS = max(MeanTS, m_MaxTS);
	uint IdxA = DA.m_ChainA->m_Idx;
	uint IdxB = DA.m_ChainB->m_Idx;
	asserta(IdxA < m_ChainCount);
	asserta(IdxB < m_ChainCount);
	fprintf(m_fDistMx, "%u\t%u\t%.3f\n", IdxA, IdxB, MeanTS);
	}

void cmd_distmx()
	{
	asserta(optset_output);
	const string &DBFN = g_Arg1;

	DistMxSearcher DBS;
	DSSParams Params;
	Params.SetFromCmdLine(10000);
	DBS.m_Params = &Params;
	DBS.LoadDB(DBFN);

	Params.m_DBSize = (float) DBS.GetDBSize();
	DBS.m_fDistMx = CreateStdioFile(opt_output);
	DBS.Setup();
	DBS.RunSelf();
	ProgressLog("maxts %.3f\n", DBS.m_MaxTS);
	CloseStdioFile(DBS.m_fDistMx);
	}
