#pragma once

#include "alpha.h"
#include "kappa_dex.h"
#include "kappa_mermx.h"
#include "mudex.h"
#include "twohitdiag.h"
#include "diaghsp.h"
#include "seqdb.h"
#include "diag.h"
#include "rankedscoresbag.h"
#include "kappa_filter_params.h"
#include "kappa_seqsource.h"

extern int16_t kappa32_flat_logodds[32*32];

#define	TRACE			0

///////////////////////////////////////////////////////////
// For one target sequence build a list of query sequences
// with 2-kmer diagonals and their scores.
///////////////////////////////////////////////////////////

class kappa_filter
	{
public:
	static RankedScoresBag m_RSB;
	// Batched RSB updates: flush when this many entries are pending.
	static const uint RSB_BATCH = 512;

	static uint8_t **m_query_kappa_codeseq_vec;
	static const uint *m_query_lengths;
	//static uint m_NQ;
	static kappa_seqsource *m_db_seqsource;
	static atomic<time_t> m_time_last_progress;
	static const kappa_mermx *m_ptrScoreMx;
	static const kappa_dex *m_ptrQKmerIndex;

#if TRACE
public:
	uint m_Trace_QIdx = UINT_MAX;
	uint m_Trace_TIdx = UINT_MAX;
	uint m_TBaseKmer = UINT_MAX;
	bool DoTrace(uint QIdx) const
		{ 
		return true;
		}
#endif

public:
///////////////////////////////////////////////////
// Query DB is typically smaller, indexed in memory
// Sequences are integers 0..19 not ASCII chars
///////////////////////////////////////////////////
	static uint m_QSeqCount;

//////////////////////////////////////
// Index of k-mers in the Query.
//////////////////////////////////////
	const kappa_dex *m_QKmerIndex = 0;
	const int16_t *m_KmerSelfScores = 0;

	const kappa_mermx *m_ScoreMx = 0;

///////////////////////////////////////////////////
// Accumulating results for current target sequence
///////////////////////////////////////////////////
	uint m_NrQueriesWithTwoHitDiag = UINT_MAX;
	uint16_t *m_QSeqIdxsWithTwoHitDiag = 0;
	uint16_t *m_QSeqIdxToBestDiagScore = 0;

// High-scoring k-mers in the neighborhood
// of the current Target k-mer
	uint *m_NeighborKmers = 0;

//////////////////////////////////////////////////////
//  m_DiagBag stores k-mer matches between the current
//  Target sequence and Query sequences.
// 	Matches are stored as (QSeqIdx, DiagIdx) pairs.
//  After scanning all k-mers, TwoHitDiag::SetUniqueFine()
//  dedupes fine (QSeqIdx, Diag) pairs before HSP extension.
//////////////////////////////////////////////////////
	TwoHitDiag m_DiagBag;

//////////////////////////////////////////////////////
// Current Target sequence
//////////////////////////////////////////////////////
	uint m_TSeqIdx = UINT_MAX;
	const byte *m_TSeq = 0;
	string m_TLabel;
	uint m_TL = UINT_MAX;
	vector<uint> m_TKmers;

	// Pending (query, target, score) for batched 
	//   AddScoresBatch; not cleared per target.
	vector<RankedScoreBatchEntry> m_RSBPending;

public:
	kappa_filter() = default;
	~kappa_filter();

	void SetQDB(const SeqDB &QDB);
	void Search(uint TSeqIdx, const string &TLabel,
				const byte *TSeq, uint TL);
	void Search_TargetSeq(uint TSeqIdx, const string &TLabel,
				   const byte *TSeq, uint TL);
	int FindHSP(uint QSeqIdx, int Diag) const;
	int FindHSP(const byte *QSeq, uint QL, int Diag) const;
	int FindHSP2(uint QSeqIdx, int Diag, int &Lo, int &Len) const;
	int FindHSP2(const byte *QSeq, uint QL, int Diag, int &Lo, int &Len) const;
	void Search_TargetKmers();
	void Search_TargetKmerNeighborhood(uint Kmer, uint TPos);
	void Search_TargetKmer(uint Kmer, uint TPos);
	void FindTwoHitDiags();
	void ExtendTwoHitDiagsToHSPs();
	int ExtendDiagToHSP(uint32_t QSeqIdx, uint16_t Diag);
	void AddTwoHitDiag(uint QSeqIdx, uint16_t Diag, int DiagScore);
	void GetResults(vector<uint> &QSeqIdxs,
					vector<uint16_t> &DiagScores) const;
	void LogDiag(uint QSeqIdx, uint16_t Diag) const;
	const char *KmerToStr(uint Kmer, string &s) const
		{
		return m_QKmerIndex->KmerToStr(Kmer, s);
		}
	void Reset();
	void LogQueryKmers(uint QSeqIdx) const;
	uint GetQKmer(uint QSeqIdx, uint QPos) const;
	void LogTargetKmers() const;
	void alloc();

public:
	static bool m_init_kappa_done;

public:
	static void init_kappa();
	static void run_filter(
		uint8_t **query_kappa_codeseqs,
		const uint *query_lengths,
		uint NQ, kappa_seqsource &db_ss);
	static void static_thread_body(uint threadidx);
	static void static_bcb_thread_body(uint threadidx);
	};
