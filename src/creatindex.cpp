#include "myutils.h"
#include "kappa_dex.h"
#include "kappa_seqsource.h"
#include "kappa_mermx.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "objmgr.h"
#include "seqinfo.h"

static void RoundTripCheck(const kappa_dex &A, const kappa_dex &B)
	{
	asserta(A.m_k == B.m_k);
	asserta(A.m_K == B.m_K);
	asserta(A.m_DictSize == B.m_DictSize);
	asserta(A.m_nseq == B.m_nseq);
	asserta(A.m_Size == B.m_Size);
	asserta(A.m_MinKmerSelfScore == B.m_MinKmerSelfScore);
	asserta(A.m_Offsets != 0 && B.m_Offsets != 0);
	for (uint i = 0; i < A.m_k; ++i)
		asserta(A.m_Offsets[i] == B.m_Offsets[i]);
	asserta(A.m_Finger != 0 && B.m_Finger != 0);
	asserta(A.m_RowSizes != 0 && B.m_RowSizes != 0);
	for (uint i = 0; i < A.m_DictSize + 2; ++i)
		asserta(A.m_Finger[i] == B.m_Finger[i]);
	for (uint i = 0; i < A.m_DictSize; ++i)
		asserta(A.m_RowSizes[i] == B.m_RowSizes[i]);
	const uint64 DataBytes = uint64(A.m_Size) * uint64(kappa_dex::m_ItemSize);
	if (DataBytes > 0)
		{
		asserta(A.m_Data != 0 && B.m_Data != 0);
		asserta(memcmp(A.m_Data, B.m_Data, size_t(DataBytes)) == 0);
		}
	ProgressLog("Round-trip OK\n");
	}

static void OpenKSS(kappa_seqsource &KSS, BCAData &bcb)
	{
	if (EndsWith(g_Arg1, ".bcb"))
		{
		bcb.Open(g_Arg1);
		KSS.OpenBCB(bcb);
		}
	else
		KSS.OpenFasta(g_Arg1);

	}

void cmd_createindex()
	{
	asserta(optset_output);
	kappa_seqsource KSS;
	BCAData bcb;
	OpenKSS(KSS, bcb);

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
	ProgressStep(999, 1000, "Pass 1");

	KmerIndex.AdjustFinger();
	KmerIndex.Alloc_Pass2();
	idx = 0;
	bcb.Close();
	KSS.Close();

	OpenKSS(KSS, bcb);
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
	ProgressStep(999, 1000, "Pass 2");
	KmerIndex.m_nseq = idx;
	KmerIndex.SetRowSizes();
	KmerIndex.LogStats();

	KmerIndex.ToFile(opt(output));

	kappa_dex Check;
	Check.FromFile(opt(output));
	RoundTripCheck(KmerIndex, Check);
	Check.LogStats();
	}

void cmd_idx_stats()
	{
	kappa_dex Index;
	Index.FromFile(g_Arg1);
	ProgressLog("k=%u K=%u dict=%u nseq=%u postings=%u minself=%d\n",
		Index.m_k, Index.m_K, Index.m_DictSize, Index.m_nseq, Index.m_Size,
		Index.m_MinKmerSelfScore);
	Index.LogStats();
	}
