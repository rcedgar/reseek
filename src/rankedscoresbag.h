#pragma once

///////////////////////////////////////
// RankedScoresBag is a container for
// sorted lists of high-scoring
// Target indexes of maximum size B.
// If there are >B targets, then only
// the top B are kept.
// There is one list per Query.
///////////////////////////////////////

#define	CHECK_SCORE_VECS	0
#define STORE_PAIR_SCORES	1	// TODO

struct RankedScoreBatchEntry
	{
	uint QueryIdx;
	uint TargetIdx;
	uint16_t Score;
	};

class RankedScoresBag
	{
public:
	uint m_B = 0;
	vector<vector<uint16_t> > m_QueryIdxToScoreVec;
	vector<vector<uint> > m_QueryIdxToTargetIdxVec;
	vector<uint16_t> m_QueryIdxToLoScore;
	uint m_QueryCount = UINT_MAX;

	mutex m_DataLock;
#if CHECK_SCORE_VECS
	vector<vector<uint16_t> > m_QueryIdxToFullScoreVec;
	vector<vector<uint> > m_QueryIdxToFullTargetIdxVec;
#endif
	vector<vector<uint16_t> > m_QueryIdxToTopScoreVec;

public:
	void Init(uint QueryCount);
	uint TruncateAllQueryVecs();
	void TruncateVecs(uint QIdx);
	const vector<uint> &GetTargetIdxs(uint QueryIdx) const;
	void GetTargetInfo(vector<uint> &TargetIdxs,
		unordered_map<uint, vector<uint> > &TargetIdxToQueryIdxs) const;
	void GetTargetInfoSorted(vector<uint> &TargetIdxs,
		unordered_map<uint, vector<uint> > &TargetIdxToQueryIdxs,
		unordered_map<uint, vector<uint> > &TargetIdxToDiagScores,
		uint &max_queries_per_target) const;

	void AddScore(uint QueryIdx, uint TargetIdx, uint16_t Score);
	// Caller must hold m_DataLock (e.g. use AddScoresBatch for batched updates).
	void AddScore_unlocked(uint QueryIdx, uint TargetIdx, uint16_t Score);
	// Sorts by QueryIdx, applies under one lock, clears Batch.
	void AddScoresBatch(vector<RankedScoreBatchEntry> &Batch);
	void ToTsv(FILE *fTsv);
	void ToLabelsTsv(FILE *fTsv,
					 const vector<string> &QLabels,
					 const vector<string> &TLabels);
#if CHECK_SCORE_VECS
	void CheckScoreVecs(uint QIdx);
	void CheckAllScoreVecs();
#endif
	};
