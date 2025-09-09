#pragma once

#include "alpha.h"
#include "mermx.h"
#include "mudex.h"
#include "twohitdiag.h"
#include "diaghsp.h"
#include "seqdb.h"
#include "diag.h"
#include "rankedscoresbag.h"
#include "prefiltermuparams.h"

const MerMx &GetMuMerMx(uint k);

extern int8_t Mu_S_ij_i8[36][36];

#define	TRACE			0

///////////////////////////////////////////////////////////
// For one target sequence build a list of query sequences
// with 2-kmer diagonals and their scores.
///////////////////////////////////////////////////////////

class PrefilterMu
	{
public:
	static mutex m_Lock;
	static RankedScoresBag m_RSB;

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
	const SeqDB *m_QDB = 0;
	uint m_QSeqCount = UINT_MAX;

//////////////////////////////////////
// Index of k-mers in the Query.
// These are 3Di 6-mers
//////////////////////////////////////
	const MuDex *m_QKmerIndex = 0;
	const int16_t *m_KmerSelfScores = 0;

////////////////////////////////////////////
//  ScoreMx is the Mu substitution matrix
//  with support for enumerating k-mers in
//  the high-scoring neighborhood of a k-mer
////////////////////////////////////////////
	const MerMx *m_ScoreMx = 0;

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
//  After scanning all k-mers, TwoHitDiag::SetDupes()
//  finds two-hit diagonals.
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
#if KMER_SORT
	vector<uint> m_TKmerSizes;
	vector<uint> m_TKmerSizeOrder;
#endif

public:
	void SetQDB(const SeqDB &QDB);
	void Search(FILE *fTsv, uint TSeqIdx, const string &TLabel,
				const byte *TSeq, uint TL);
	void SetTarget(uint TSeqIdx, const string &TLabel,
				   const byte *TSeq, uint TL);
	void Search_TargetSeq();
	int FindHSP(uint QSeqIdx, int Diag) const;
	int FindHSP2(uint QSeqIdx, int Diag, int &Lo, int &Len) const;
	void Search_TargetKmers();
	void Search_TargetKmerNeighborhood(uint Kmer, uint TPos);
	void Search_Kmer(uint Kmer, uint TPos);
	void FindTwoHitDiags();
	void ExtendTwoHitDiagsToHSPs();
	int ExtendTwoHitDiagToHSP(uint32_t QSeqIdx, uint16_t Diag);
	void AddTwoHitDiag(uint QSeqIdx, uint16_t Diag, int DiagScore);
	void GetResults(vector<uint> &QSeqIdxs,
					vector<uint16_t> &DiagScores) const;
	void LogDiag(uint QSeqIdx, uint16_t Diag) const;
	const char *KmerToStr(uint Kmer, string &s) const
		{
		return m_QKmerIndex->KmerToStr(Kmer, s);
		}
	void Reset();
	};
