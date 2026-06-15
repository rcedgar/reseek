#include "myutils.h"
#include "rankedscoresbag.h"
#include "sort.h"
#include <algorithm>

const vector<uint> &RankedScoresBag::GetTargetIdxs(uint QueryIdx) const
	{
	asserta(QueryIdx < SIZE(m_QueryIdxToTargetIdxVec));
	return m_QueryIdxToTargetIdxVec[QueryIdx];
	}

void RankedScoresBag::TruncateVecs(uint QueryIdx)
	{
	vector<uint16_t> &ScoreVec = m_QueryIdxToScoreVec[QueryIdx];
	uint CurrentSize = SIZE(ScoreVec);
	if (CurrentSize < m_B)
		return;

	uint *Order = myalloc(uint, CurrentSize);
	QuickSortOrderDesc(ScoreVec.data(), CurrentSize, Order);
	vector<uint> &TargetIdxVec = m_QueryIdxToTargetIdxVec[QueryIdx];
	vector<uint16_t> NewScoreVec;
	vector<uint> NewTargetIdxVec;
	NewScoreVec.reserve(m_B);
	NewTargetIdxVec.reserve(m_B);
	for (uint k = 0; k < m_B; ++k)
		{
		uint i = Order[k];
		uint16_t Score = ScoreVec[i];
		uint TargetIdx = TargetIdxVec[i];
		NewScoreVec.push_back(Score);
		NewTargetIdxVec.push_back(TargetIdx);
		}
	uint16_t NewLo = NewScoreVec[m_B-1];
	m_QueryIdxToTargetIdxVec[QueryIdx] = NewTargetIdxVec;
	m_QueryIdxToScoreVec[QueryIdx] = NewScoreVec;
	m_QueryIdxToLoScore[QueryIdx] = NewLo;
	myfree(Order);
	}

void RankedScoresBag::AddScore_unlocked(uint QueryIdx, uint TargetIdx, uint16_t Score)
	{
#if CHECK_SCORE_VECS
	m_QueryIdxToFullTargetIdxVec[QueryIdx].push_back(TargetIdx);
	m_QueryIdxToFullScoreVec[QueryIdx].push_back(Score);
#endif
	vector<uint16_t> &ScoreVec = m_QueryIdxToScoreVec[QueryIdx];
	uint16_t LoScore = m_QueryIdxToLoScore[QueryIdx];
	if (Score >= LoScore)
		{
		ScoreVec.push_back(Score);
		m_QueryIdxToTargetIdxVec[QueryIdx].push_back(TargetIdx);
		if (SIZE(ScoreVec) >= 2*m_B)
			TruncateVecs(QueryIdx);
#if STORE_PAIR_SCORES
		{
		asserta(QueryIdx < m_QueryIdxToTopScoreVec.size());
		asserta(TargetIdx < m_QueryIdxToTopScoreVec[QueryIdx].size());
		uint16_t TopScore =
			m_QueryIdxToTopScoreVec[QueryIdx][TargetIdx];
		if (Score > TopScore)
			m_QueryIdxToTopScoreVec[QueryIdx][TargetIdx] = Score;
		}
#endif
		}
	}

void RankedScoresBag::AddScore(uint QueryIdx, uint TargetIdx, uint16_t Score)
	{
	m_DataLock.lock();
	AddScore_unlocked(QueryIdx, TargetIdx, Score);
	m_DataLock.unlock();
	}

void RankedScoresBag::AddScoresBatch(vector<RankedScoreBatchEntry> &Batch)
	{
	const uint N = SIZE(Batch);
	if (N == 0)
		return;
	std::sort(Batch.begin(), Batch.end(),
		 [](const RankedScoreBatchEntry &a, const RankedScoreBatchEntry &b)
			{
			return a.QueryIdx < b.QueryIdx;
			});
	m_DataLock.lock();
	for (uint i = 0; i < N; ++i)
		{
		const RankedScoreBatchEntry &e = Batch[i];
		AddScore_unlocked(e.QueryIdx, e.TargetIdx, e.Score);
		}
	m_DataLock.unlock();
	Batch.clear();
	}

#if CHECK_SCORE_VECS
void RankedScoresBag::CheckScoreVecs(uint QueryIdx)
	{
	TruncateVecs(QueryIdx);

	const vector<uint16_t> &ScoreVec = m_QueryIdxToScoreVec[QueryIdx];
	const vector<uint> &TargetIdxVec = m_QueryIdxToTargetIdxVec[QueryIdx];
	const uint N = SIZE(ScoreVec);
	if (N == 0)
		{
		asserta(ScoreVec.empty());
		asserta(TargetIdxVec.empty());
		return;
		}
	uint16_t LoScore = m_QueryIdxToLoScore[QueryIdx];
	asserta(SIZE(TargetIdxVec) == N);
	map<uint, uint16_t> TargetIdxToScore;
	for (uint i = 0; i < N; ++i)
		{
		uint16_t Score = ScoreVec[i];
		asserta(Score >= LoScore);
		uint TargetIdx = TargetIdxVec[i];
		asserta(TargetIdxToScore.find(TargetIdx) == TargetIdxToScore.end());
		TargetIdxToScore[TargetIdx] = Score;
		}

	const vector<uint16_t> &FullScoreVec = m_QueryIdxToFullScoreVec[QueryIdx];
	const vector<uint> &FullTargetIdxVec = m_QueryIdxToFullTargetIdxVec[QueryIdx];
	uint FullSize = SIZE(FullScoreVec);
	asserta(SIZE(FullTargetIdxVec) == FullSize);
	uint *Order = myalloc(uint, FullSize);
	QuickSortOrderDesc(FullScoreVec.data(), FullSize, Order);
	uint M = FullSize - 1;
	if (M > m_B-1)
		M = m_B-1;
	uint Lok = Order[M];
	int16_t FullLoScore = FullScoreVec[Lok];

	for (uint k = 0; k < FullSize; ++k)
		{
		uint i = Order[k];
		assert(i < FullSize);
		uint TargetIdx = FullTargetIdxVec[i];
		uint16_t Score = FullScoreVec[i];
		map<uint, uint16_t>::const_iterator iter = TargetIdxToScore.find(TargetIdx);
		if (Score > LoScore)
			{
			if (iter == TargetIdxToScore.end())
				{
				Log("CheckScoreVecs(%u) LoScore=%d, FullLoScore=%d\n",
					QueryIdx, LoScore, FullLoScore);
				Log("Short: ");
				for (uint i = 0; i < N; ++i)
					Log(" %u=%d", TargetIdxVec[i], ScoreVec[i]);
				Log("\n");
				Log("Full: ");
				for (uint k = 0; k < FullSize; ++k)
					{
					uint i = Order[k];
					Log(" %u=%d", FullTargetIdxVec[i], FullScoreVec[i]);
					}
				Log("\n");
				Die("CheckAllScoreVecs");
				}
			asserta(iter->second == Score);
			}
		}
	myfree(Order);
	}

void RankedScoresBag::CheckAllScoreVecs()
	{
	const uint QueryCount = SIZE(m_QueryIdxToScoreVec);
	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		CheckScoreVecs(QueryIdx);
	}
#endif // CHECK_SCORE_VECS

void RankedScoresBag::ToLabelsTsv(FILE *f,
					const vector<string> &QLabels,
					const vector<string> &TLabels)
	{
	if (f == 0)
		return;

	asserta(SIZE(QLabels) == m_QueryCount);
	const uint TCount = SIZE(TLabels);

	for (uint QueryIdx = 0; QueryIdx < m_QueryCount; ++QueryIdx)
		TruncateVecs(QueryIdx);

	map<uint, vector<uint> > TargetIdxToQueryIdxs;
	vector<uint> TargetIdxs;
	for (uint QueryIdx = 0; QueryIdx < m_QueryCount; ++QueryIdx)
		{
		const vector<uint16_t> &ScoreVec = m_QueryIdxToScoreVec[QueryIdx];
		const vector<uint> &TargetIdxVec = m_QueryIdxToTargetIdxVec[QueryIdx];
		const uint n = SIZE(ScoreVec);
		for (uint i = 0; i < n; ++i)
			{
			uint TargetIdx = TargetIdxVec[i];
			if (TargetIdxToQueryIdxs.find(TargetIdx) == TargetIdxToQueryIdxs.end())
				{
				TargetIdxs.push_back(TargetIdx);
				vector<uint> Empty;
				TargetIdxToQueryIdxs[TargetIdx] = Empty;
				}
			TargetIdxToQueryIdxs[TargetIdx].push_back(QueryIdx);
			}
		}
	const uint TargetCount = SIZE(TargetIdxs);
	QuickSortInPlace(TargetIdxs.data(), TargetCount);
	for (uint k = 0; k < TargetCount; ++k)
		{
		uint TargetIdx = TargetIdxs[k];
		asserta(TargetIdx < TCount);
		const string &TLabel = TLabels[TargetIdx];

		map<uint, vector<uint> >::const_iterator iter = TargetIdxToQueryIdxs.find(TargetIdx);
		asserta(iter != TargetIdxToQueryIdxs.end());
		const vector<uint> &QIdxs = iter->second;
		const uint K = SIZE(QIdxs);
		for (uint i = 0; i < K; ++i)
			{
			uint QIdx = QIdxs[i];
			asserta(QIdx < m_QueryCount);
			const string &QLabel = QLabels[QIdx];
			fprintf(f, "%s\t%s\n", QLabel.c_str(), TLabel.c_str());
			}
		}
	}

uint RankedScoresBag::TruncateAllQueryVecs()
	{
	uint Total = 0;
	for (uint QueryIdx = 0; QueryIdx < m_QueryCount; ++QueryIdx)
		{
		ProgressStep(QueryIdx, m_QueryCount, "Truncate prefilter vecs");
		TruncateVecs(QueryIdx);
		Total += SIZE(m_QueryIdxToTargetIdxVec[QueryIdx]);
		}
	return Total;
	}

void RankedScoresBag::GetTargetInfo(
	vector<uint> &TargetIdxs,
	unordered_map<uint, vector<uint> > &TargetIdxToQueryIdxs) const
	{
	TargetIdxToQueryIdxs.clear();
	TargetIdxs.clear();
	Progress("Make target info (unsorted)...\n");
	for (uint QueryIdx = 0; QueryIdx < m_QueryCount; ++QueryIdx)
		{
		const vector<uint16_t> &ScoreVec = m_QueryIdxToScoreVec[QueryIdx];
		const vector<uint> &TargetIdxVec = m_QueryIdxToTargetIdxVec[QueryIdx];
		const uint n = SIZE(ScoreVec);
		for (uint i = 0; i < n; ++i)
			{
			uint TargetIdx = TargetIdxVec[i];
			if (TargetIdxToQueryIdxs.find(TargetIdx) == TargetIdxToQueryIdxs.end())
				{
				TargetIdxs.push_back(TargetIdx);
				vector<uint> Empty;
				TargetIdxToQueryIdxs[TargetIdx] = Empty;
				}
			TargetIdxToQueryIdxs[TargetIdx].push_back(QueryIdx);
			}
		}
	const uint TargetCount = SIZE(TargetIdxs);
	QuickSortInPlace(TargetIdxs.data(), TargetCount);
	Progress("%u targets\n", TargetCount);
	}

void RankedScoresBag::GetTargetInfoSorted(
	vector<uint> &TargetIdxs,
	unordered_map<uint, vector<uint> > &TargetIdxToQueryIdxs,
	unordered_map<uint, vector<uint> > &TargetIdxToDiagScores,
	uint &max_queries_per_target) const
	{
	max_queries_per_target = 0;
	TargetIdxToQueryIdxs.clear();
	TargetIdxToDiagScores.clear();
	TargetIdxs.clear();
	Progress("Make target info (sorted)...\n");
	unordered_map<uint, vector<pair<uint, uint16_t> > > TargetIdxToPairs;
	for (uint QueryIdx = 0; QueryIdx < m_QueryCount; ++QueryIdx)
		{
		const vector<uint16_t> &ScoreVec = m_QueryIdxToScoreVec[QueryIdx];
		const vector<uint> &TargetIdxVec = m_QueryIdxToTargetIdxVec[QueryIdx];
		const uint n = SIZE(ScoreVec);
		for (uint i = 0; i < n; ++i)
			{
			uint TargetIdx = TargetIdxVec[i];
			uint16_t Score = ScoreVec[i];
			if (TargetIdxToPairs.find(TargetIdx) == TargetIdxToPairs.end())
				{
				TargetIdxs.push_back(TargetIdx);
				vector<pair<uint, uint16_t> > Empty;
				TargetIdxToPairs[TargetIdx] = Empty;
				}
			TargetIdxToPairs[TargetIdx].push_back(make_pair(QueryIdx, Score));
			}
		}
	const uint TargetCount = SIZE(TargetIdxs);
	QuickSortInPlace(TargetIdxs.data(), TargetCount);
	for (uint k = 0; k < TargetCount; ++k)
		{
		uint TargetIdx = TargetIdxs[k];
		vector<pair<uint, uint16_t> > &Pairs = TargetIdxToPairs[TargetIdx];
		std::sort(Pairs.begin(), Pairs.end(),
			[](const pair<uint, uint16_t> &a, const pair<uint, uint16_t> &b)
				{
				return a.second > b.second;
				});
		vector<uint> &QIdxs = TargetIdxToQueryIdxs[TargetIdx];
		vector<uint> &DiagScores = TargetIdxToDiagScores[TargetIdx];
		uint nq = uint(Pairs.size());
		max_queries_per_target = max(nq, max_queries_per_target);
		QIdxs.reserve(nq);
		DiagScores.reserve(nq);
		for (uint i = 0; i < nq; ++i)
			{
			QIdxs.push_back(Pairs[i].first);
			DiagScores.push_back(Pairs[i].second);
			}
		}
	Progress("%u targets, maxqpert %u\n", TargetCount, max_queries_per_target);
	}

void RankedScoresBag::ToTsv(FILE *f)
	{
	if (f == 0)
		return;

	unordered_map<uint, vector<uint> > TargetIdxToQueryIdxs;
	vector<uint> TargetIdxs;
	GetTargetInfo(TargetIdxs, TargetIdxToQueryIdxs);

	const uint TargetCount = SIZE(TargetIdxs);

	fprintf(f, "prefilter\t%u\n", TargetCount);
	uint64_t Total = 0;
	for (uint k = 0; k < TargetCount; ++k)
		{
		uint TargetIdx = TargetIdxs[k];
		unordered_map<uint, vector<uint> >::const_iterator iter =
			TargetIdxToQueryIdxs.find(TargetIdx);
		asserta(iter != TargetIdxToQueryIdxs.end());
		const vector<uint> &QIdxs = iter->second;
		const uint K = SIZE(QIdxs);
		Total += K;
		fprintf(f, "%u\t%u", TargetIdx, K);
		for (uint i = 0; i < K; ++i)
			fprintf(f, "\t%u", QIdxs[i]);
		fprintf(f, "\n");
		}
	ProgressLog("%s prefilter hits\n", FloatToStr(double(Total)));
	}

void RankedScoresBag::Init(uint QueryCount)
	{
	m_QueryCount = QueryCount;
	m_QueryIdxToScoreVec.resize(QueryCount);
	m_QueryIdxToTargetIdxVec.resize(QueryCount);
	m_QueryIdxToLoScore.resize(QueryCount, 0);
#if CHECK_SCORE_VECS
	m_QueryIdxToFullTargetIdxVec.resize(QueryCount);
	m_QueryIdxToFullScoreVec.resize(QueryCount);
#endif
#if STORE_PAIR_SCORES
	m_QueryIdxToTopScoreVec.resize(QueryCount);
	for (uint i = 0; i < QueryCount; ++i)
		m_QueryIdxToTopScoreVec[i].resize(QueryCount);
#endif
	}
