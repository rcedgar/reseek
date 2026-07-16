#include "myutils.h"
#include "kappa_dex.h"
#include "kappa_seqsource.h"
#include "kappa_mermx.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "objmgr.h"
#include "seqinfo.h"

void cmd_createindex()
	{
	kappa_seqsource KSS;
	KSS.OpenFasta(g_Arg1);

	flat_params params;
	params.init_from_cmdline();
	params.logme();

	kappa_dex KmerIndex;
	KmerIndex.Init();

	const uint k = flat_params::m_kappa_kmer_nrones;
	const kappa_mermx &GetKappaMerMx(uint k);
	const kappa_mermx &ScoreMx = GetKappaMerMx(k);
	asserta(ScoreMx.m_k == k);

	KmerIndex.m_KmerSelfScores = ScoreMx.BuildSelfScores_Kmers();
	KmerIndex.m_MinKmerSelfScore = flat_params::m_kappa_min_kmerpairscore;
	KmerIndex.m_AddNeighborhood = false;
	KmerIndex.m_ptrScoreMx = 0;

	KmerIndex.Alloc_Pass1();
	
	ObjMgr OM;
	SeqInfo *SI = OM.GetSeqInfo();
	ProgressStep(0, 1000, "Pass 1");
	uint idx = 0;
	for (;;)
		{
		uint PctX10 = KSS.GetPctDoneX10();
		bool ok = KSS.GetNext(SI);
		if (!ok)
			break;
		if (PctX10 > 0) ProgressStep(PctX10, 1000, "Pass 1");
		KmerIndex.SetSeq(idx++, SI->m_Label, SI->m_Seq, SI->m_L);
		KmerIndex.AddSeq_Pass1();
		}
	ProgressStep(0, 1000, "Pass 1");

	KmerIndex.AdjustFinger();
	KmerIndex.Alloc_Pass2();
	idx = 0;
	KSS.Close();
	KSS.OpenFasta(g_Arg1);
	for (;;)
		{
		uint PctX10 = KSS.GetPctDoneX10();
		bool ok = KSS.GetNext(SI);
		if (!ok)
			break;
		if (PctX10 > 0) ProgressStep(PctX10, 1000, "Pass 2");
		KmerIndex.SetSeq(idx++, SI->m_Label, SI->m_Seq, SI->m_L);
		KmerIndex.AddSeq_Pass2();
		}
	ProgressStep(0, 1000, "Pass 2");
	KmerIndex.SetRowSizes();
	KmerIndex.LogStats();
	}
