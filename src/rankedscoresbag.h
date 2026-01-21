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

public:
	void Init(uint QueryCount);
	void TruncateVecs(uint QIdx);
	void AddScore(uint QueryIdx, uint TargetIdx, uint16_t Score);
	void ToTsv(FILE *fTsv);
	void ToLabelsTsv(FILE *fTsv,
					 const vector<string> &QLabels,
					 const vector<string> &TLabels);
	void ToScoreTsv(FILE *fTsv,
					 const vector<string> &QLabels,
					 const vector<string> &TLabels);
#if CHECK_SCORE_VECS
	void CheckScoreVecs(uint QIdx);
	void CheckAllScoreVecs();
#endif
	};
