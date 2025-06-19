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
	if (DA.m_EvalueA > m_MaxEvalue)
		return;
	asserta(m_fDistMx!= 0);
	uint IdxA = DA.m_ChainA->m_Idx;
	uint IdxB = DA.m_ChainB->m_Idx;
	asserta(IdxA < m_ChainCount);
	asserta(IdxB < m_ChainCount);
	float ts = DA.m_NewTestStatisticA;
	m_MaxTS = max(ts, m_MaxTS);
	fprintf(m_fDistMx, "%u\t%u\t%.3f\n", IdxA, IdxB, ts);
	}

void cmd_distmx()
	{
	asserta(optset_output);
	const string &DBFN = g_Arg1;

	if (!optset_fast && !optset_sensitive && !optset_verysensitive)
		{
		opt(fast) = true;
		optset_fast = true;
		}

	DistMxSearcher DBS;
	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast);
	DBS.m_Params = &Params;
	DBS.LoadDB(DBFN);

	DBS.m_fDistMx = CreateStdioFile(opt(output));
	DBS.Setup();
	DBS.RunSelf();
	ProgressLog("maxts %.3f\n", DBS.m_MaxTS);
	CloseStdioFile(DBS.m_fDistMx);
	}
