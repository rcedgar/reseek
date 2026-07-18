#include "myutils.h"
#include "kappa_dex.h"
#include "kappa_seqsource.h"
#include "kappa_mermx.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "sort.h"
#include "bcadata.h"

//// Debug: -bail N → exit(0) at stage N (skips later code and local dtors;
//// still runs CRT atexit / heap teardown). Stages printed in ProgressLog.
//static void BailAt(unsigned Stage, const char *Msg)
//	{
//	if (!optset_bail || opt_bail != Stage)
//		return;
//	ProgressLog("BAIL %u: %s\n", Stage, Msg);
//	if (g_fLog != 0)
//		fflush(g_fLog);
//	fflush(stdout);
//	fflush(stderr);
//	exit(0);
//	}

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

static void BuildIndexFromFasta(kappa_dex &KmerIndex, const string &FN)
	{
	kappa_seqsource KSS;
	KSS.OpenFasta(FN);

	KmerIndex.Alloc_Pass1();
	ObjMgr OM;
	SeqInfo *SI = OM.GetSeqInfo();
	uint idx = 0;
	ProgressStep(0, 1000, "Pass 1");
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
	KSS.Close();
	KSS.OpenFasta(FN);
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
	KSS.Close();
	}

void set_default_stats()
	{
	if (!optset_stats) opt_stats = mystrsave("sf");
	optset_stats = true;
	optset_fast = true;
	opt_fast = true;
	}

void cmd_createindex()
	{
	asserta(optset_output);
	set_default_stats();

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

	if (EndsWith(g_Arg1, ".bcb"))
		{
		BCAData bcb;
		bcb.Open(g_Arg1);
		asserta(bcb.m_HasNuSequences);

		uint8_t **kappa_codeseqs = 0;
		uint *lengths = 0;
		bcb.make_kappa_codeseqs(&kappa_codeseqs, &lengths);

		const uint nseq = bcb.GetChainCount();
		KmerIndex.from_codeseqs(kappa_codeseqs, lengths, bcb.m_Labels, nseq);

		bcb.Close();
		}
	else
		{
		BuildIndexFromFasta(KmerIndex, g_Arg1);
		}

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
	ProgressLog("k=%u K=%u dict=%u nseq=%u postings=%s minself=%d\n",
		Index.m_k, Index.m_K, Index.m_DictSize, Index.m_nseq,
		Int64ToStr(Index.m_Size),
		Index.m_MinKmerSelfScore);
	Index.LogStats();
	}
