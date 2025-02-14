#pragma once

#include "alpha.h"
#include "mermx.h"
#include "threedex.h"
#include "twohitdiag.h"
#include "diaghsp.h"
#include "seqdb.h"
#include "diag.h"

///////////////////////////////////////////////////////////
// For each query sequence, build a list of at most maxseqs
// Target sequences sorted by decreasing diagonal score.
///////////////////////////////////////////////////////////

class Prefilter
	{
public:
	int m_KmerNeighborMinScore = 78;

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
	const ThreeDex *m_QKmerIndex = 0;

////////////////////////////////////////////
//  ScoreMx is the 3Di substitution matrix
//  with support for enumerating k-mers in
//  the high-scoring neighborhood of a k-mer
////////////////////////////////////////////
	const MerMx *m_ScoreMx = 0;

///////////////////////////////////////////////////
// Accumulating results for current target sequence
///////////////////////////////////////////////////
	uint m_NrQueriesWithTwoHitDiag = UINT_MAX;
	uint *m_QSeqIdxsWithTwoHitDiag = 0;
	int *m_QSeqIdxToBestDiagScore = 0;

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

	const byte *m_TSeq = 0;			// current Target sequence
	uint m_TL = UINT_MAX;			// length of m_TSeq
	uint m_TSeqPos = UINT_MAX;		// position of Target k-mer

public:
	void SetQDB(const SeqDB &QDB);
	void Search_TargetSeq(const byte *TSeq, uint TL,
						  vector<uint> &QSeqIdxs,
						  vector<int> &DiagScores);
	int FindHSP(const byte *Q, uint QL, int Diag) const;
	void Search_TargetKmers();
	void Search_TargetKmerNeighborhood(uint Kmer);
	void Search_Kmer(uint Kmer);
	void FindTwoHitDiags();
	void ExtendTwoHitDiagsToHSPs();
	int ExtendTwoHitDiagToHSP(uint32_t QSeqIdx, uint16_t Diag);
	void AddTwoHitDiag(uint QSeqIdx, uint16_t Diag, int DiagScore);
	void GetResults(vector<uint> &QSeqIdxs,
					vector<int> &DiagScores) const;
	void Reset();
	};
