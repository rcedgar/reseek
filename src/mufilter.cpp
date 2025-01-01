#include "myutils.h"
#include "seqdb.h"
#include "mukmerfilter.h"
#include "dss.h"
#include "alpha.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "fastaseqsource.h"
#include "sort.h"

#define CHECK_SCORE_VECS	1

static atomic<uint> s_PairCount;
static atomic<uint> s_BCSCount;
static uint s_TargetCount;
static mutex s_ProgressLock;
static mutex s_DataLock;
static time_t s_last_progress;

static uint MU_FILTER_KEEPN = 5;

static vector<vector<int> > s_QueryIdxToScoreVec;
static vector<vector<uint> > s_QueryIdxToTargetIdxVec;
static vector<int> s_QueryIdxToLoScore;

static void TruncateVecs(uint QueryIdx)
	{
	vector<int> &ScoreVec = s_QueryIdxToScoreVec[QueryIdx];
	uint CurrentSize = SIZE(ScoreVec);
	if (CurrentSize < MU_FILTER_KEEPN)
		return;

	uint *Order = myalloc(uint, CurrentSize);
	QuickSortOrderDesc(ScoreVec.data(), CurrentSize, Order);

	vector<uint> &TargetIdxVec = s_QueryIdxToTargetIdxVec[QueryIdx];
	vector<int> NewScoreVec;
	vector<uint> NewTargetIdxVec;
	NewScoreVec.reserve(MU_FILTER_KEEPN);
	NewTargetIdxVec.reserve(MU_FILTER_KEEPN);
	for (uint k = 0; k < MU_FILTER_KEEPN; ++k)
		{
		uint i = Order[k];
		int Score = ScoreVec[i];
		uint TargetIdx = TargetIdxVec[i];
		NewScoreVec.push_back(Score);
		NewTargetIdxVec.push_back(TargetIdx);
		}
	s_QueryIdxToTargetIdxVec[QueryIdx] = NewTargetIdxVec;
	s_QueryIdxToScoreVec[QueryIdx] = NewScoreVec;
	s_QueryIdxToLoScore[QueryIdx] = NewScoreVec[MU_FILTER_KEEPN-1];
	myfree(Order);
	}

#if CHECK_SCORE_VECS
static vector<vector<int> > s_QueryIdxToFullScoreVec;
static vector<vector<uint> > s_QueryIdxToFullTargetIdxVec;

static void CheckScoreVecs(uint QueryIdx)
	{
	TruncateVecs(QueryIdx);

	const vector<int> &ScoreVec = s_QueryIdxToScoreVec[QueryIdx];
	const vector<uint> &TargetIdxVec = s_QueryIdxToTargetIdxVec[QueryIdx];
	const uint N = SIZE(ScoreVec);
	asserta(SIZE(TargetIdxVec) == N);
	map<uint, int> TargetIdxToScore;
	for (uint i = 0; i < N; ++i)
		{
		int Score = ScoreVec[i];
		uint TargetIdx = TargetIdxVec[i];
		asserta(TargetIdxToScore.find(TargetIdx) == TargetIdxToScore.end());
		TargetIdxToScore[TargetIdx] = Score;
		}

	const vector<int> &FullScoreVec = s_QueryIdxToFullScoreVec[QueryIdx];
	const vector<uint> &FullTargetIdxVec = s_QueryIdxToFullTargetIdxVec[QueryIdx];
	uint FullSize = SIZE(FullScoreVec);
	asserta(SIZE(FullTargetIdxVec) == FullSize);
	uint *Order = myalloc(uint, FullSize);
	QuickSortOrderDesc(FullScoreVec.data(), FullSize, Order);
	uint Lok = Order[FullSize-1];
	int LoScore = FullScoreVec[Lok];

	for (uint k = 0; k < FullSize; ++k)
		{
		uint i = Order[k];
		assert(i < FullSize);
		uint TargetIdx = FullTargetIdxVec[i];
		int Score = FullScoreVec[i];
		map<uint, int>::const_iterator iter = TargetIdxToScore.find(TargetIdx);
		if (Score > LoScore)
			{
			asserta(iter != TargetIdxToScore.end());
			asserta(iter->second == Score);
			}
		}
	myfree(Order);
	}

static void CheckAllScoreVecs()
	{
	const uint QueryCount = SIZE(s_QueryIdxToScoreVec);
	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		CheckScoreVecs(QueryIdx);
	}
#endif // CHECK_SCORE_VECS

static void AddScore(uint QueryIdx, uint TargetIdx, int Score)
	{
	s_DataLock.lock();
#if CHECK_SCORE_VECS
	s_QueryIdxToFullTargetIdxVec[QueryIdx].push_back(TargetIdx);
	s_QueryIdxToFullScoreVec[QueryIdx].push_back(Score);
#endif
	vector<int> &ScoreVec = s_QueryIdxToScoreVec[QueryIdx];
	uint CurrentSize = SIZE(ScoreVec);
	if (CurrentSize < MU_FILTER_KEEPN)
		{
		ScoreVec.push_back(Score);
		s_QueryIdxToTargetIdxVec[QueryIdx].push_back(TargetIdx);
		}
	else
		{
		if (Score >= s_QueryIdxToLoScore[QueryIdx])
			{
			if (CurrentSize < 2*MU_FILTER_KEEPN)
				{
				ScoreVec.push_back(Score);
				s_QueryIdxToTargetIdxVec[QueryIdx].push_back(TargetIdx);
				}
			else
				{
				asserta(CurrentSize == 2*MU_FILTER_KEEPN);
				TruncateVecs(QueryIdx);
				}
			}
		}
	s_DataLock.unlock();
	}

static void GetMuKmers(const string &PatternStr,
					   const vector<byte> &Letters,
					   vector<uint> &Kmers)
	{
	Kmers.clear();
	const uint PatternLength = SIZE(PatternStr);
	const uint L = SIZE(Letters);
	Kmers.reserve(L);
	for (uint Pos = 0; Pos + PatternLength <= L; ++Pos)
		{
		uint Kmer = 0;
		for (uint j = 0; j < PatternLength; ++j)
			{
			if (PatternStr[j] == '1')
				{
				assert(Pos + j < SIZE(Letters));
				byte Letter = Letters[Pos + j];
				assert(Letter < 36);
				Kmer = Kmer*36 + Letter;
				}
			}
		assert(Kmer != UINT_MAX);
		Kmers.push_back(Kmer);
		}
	}

static void ThreadBody(uint ThreadIndex, 
					   const DSSParams *ptrParams,
					   FASTASeqSource *ptrFSS, 
					   vector<MuKmerFilter *> *ptrMKFs)
	{
	const DSSParams &Params = *ptrParams;
	FASTASeqSource &FSS = *ptrFSS;
	vector<MuKmerFilter *> &MKFs = *ptrMKFs;
	const string &PatternStr = Params.m_PatternStr;
	const uint QueryCount = SIZE(MKFs);

	ObjMgr OM;

	SeqInfo *SI = OM.GetSeqInfo();
	vector<byte> MuLettersT;
	vector<uint> MuKmersT;
	for (;;)
		{
		bool Ok = FSS.GetNext(SI);
		if (!Ok)
			return;
		uint TL = SI->m_L;
		const byte *T = SI->m_Seq;
		MuLettersT.clear();
		MuLettersT.reserve(TL);
		for (uint i = 0; i < TL; ++i)
			{
			byte c = T[i];
			byte Letter = g_CharToLetterMu[c];
			asserta(Letter < 36);
			MuLettersT.push_back(Letter);
			}
		GetMuKmers(PatternStr, MuLettersT, MuKmersT);
		bool DoProgress = false;
		s_ProgressLock.lock();
		++s_TargetCount;
		if (s_TargetCount%100 == 0)
			{
			time_t now = time(0);
			if (now != s_last_progress)
				{
				DoProgress = true;
				s_last_progress = now;
				}
			}
		s_ProgressLock.unlock();
		if (DoProgress)
			Progress("%u chains scanned   \r", s_TargetCount);

		for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
			{
			MuKmerFilter &MKF = *MKFs[QueryIdx];
			MKF.Align(MuLettersT, MuKmersT);
			++s_PairCount;
			int Score = MKF.m_BestChainScore;
			if (Score > 0)
				{
				++s_BCSCount;
				AddScore(QueryIdx, SI->m_Index, Score);
#if CHECK_SCORE_VECS
				CheckScoreVecs(QueryIdx);
#endif
				}
			}
		}
	}

void cmd_mufilter()
	{
	asserta(optset_db);
	const string &QueryFN = g_Arg1;
	const string &DBFN = string(opt_db);
	FASTASeqSource FSS;
	FSS.Open(DBFN);

	DSSParams Params;
	Params.SetFromCmdLine(10000);
	const string &PatternStr = Params.m_PatternStr;

	SeqDB QueryDB;
	QueryDB.FromFasta(QueryFN);
	const uint QueryCount = QueryDB.GetSeqCount();

	s_QueryIdxToScoreVec.resize(QueryCount);
	s_QueryIdxToTargetIdxVec.resize(QueryCount);
	s_QueryIdxToLoScore.resize(QueryCount, 0);
#if CHECK_SCORE_VECS
	s_QueryIdxToFullTargetIdxVec.resize(QueryCount);
	s_QueryIdxToFullScoreVec.resize(QueryCount);
#endif

	vector<MuKmerFilter *> MKFs;
	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		{
		ProgressStep(QueryIdx, QueryCount, "Indexing query");
		MuKmerFilter *MKF = new MuKmerFilter;
		MKF->SetParams(Params);
		const uint QL = QueryDB.GetSeqLength(QueryIdx);
		const string &Q = QueryDB.GetSeq(QueryIdx);
		vector<byte> *ptrMuLetters = new vector<byte>;
		vector<uint> *ptrMuKmers = new vector<uint>;
		ptrMuLetters->reserve(QL);
		for (uint i = 0; i < QL; ++i)
			{
			byte c = Q[i];
			byte Letter = g_CharToLetterMu[c];
			asserta(Letter < 36);
			ptrMuLetters->push_back(Letter);
			}
		GetMuKmers(PatternStr, *ptrMuLetters, *ptrMuKmers);
		MKF->ResetQ();
		MKF->SetQ(ptrMuLetters, ptrMuKmers);
		MKFs.push_back(MKF);
		}

	const uint ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody, ThreadIndex, &Params, &FSS, &MKFs);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	ProgressLog("%u / %u with +ve HSP score (%.1f%%)\n",
				s_BCSCount.load(), s_PairCount.load(),
				GetPct(s_BCSCount, s_PairCount));

	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		TruncateVecs(QueryIdx);

#if CHECK_SCORE_VECS
	CheckAllScoreVecs();
#endif
	}
