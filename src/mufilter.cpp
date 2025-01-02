#include "myutils.h"
#include "seqdb.h"
#include "mukmerfilter.h"
#include "dss.h"
#include "alpha.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "fastaseqsource.h"
#include "sort.h"
#include "mudex.h"

extern int8_t IntScoreMx_Mu[36][36];
static const int X = 12;
static const int MIN_HSP_SCORE = 30;
static uint MU_FILTER_KEEPN = 1000;

#define CHECK_SCORE_VECS	0
#define CHECK_SCORE_VECSX	0

static atomic<uint> s_PairCount;
static atomic<uint> s_PveCount;
static uint s_TargetCount;
static mutex s_ProgressLock;
static mutex s_DataLock;
static time_t s_last_progress;


static vector<vector<int> > s_QueryIdxToScoreVec;
static vector<vector<uint> > s_QueryIdxToTargetIdxVec;
static vector<vector<byte> > s_QueryIdxToMuLetters;
static vector<vector<uint> > s_QueryIdxToMuKmers;
static vector<int> s_QueryIdxToLoScore;
static vector<uint> s_QueryIdxToLength;
static SeqDB *s_ptrQueryDB;

static int MuXDrop(const byte *Q, int PosQ, int LQ, 
				   const byte *T, int PosT, int LT,
				   int X)
	{
	StartTimer(MuXDrop);
#if DEBUG
	int Loi = PosQ;
	int Loj = PosT;
	int Len = 0;
#endif
	int i = PosQ;
	int j = PosT;

	int FwdScore = 0;
	int BestFwdScore = 0;
#if DEBUG
	int FwdLen = 0;
#endif
	while (i < LQ && j < LT)
		{
		byte q = Q[i++];
		byte t = T[j++];
		FwdScore += IntScoreMx_Mu[q][t];
		if (FwdScore > BestFwdScore)
			{
#if DEBUG
			FwdLen = i - PosQ;
#endif
			BestFwdScore = FwdScore;
			}
		else if (FwdScore + X < BestFwdScore)
			break;
		}

	int RevScore = 0;
	int BestRevScore = 0;
#if DEBUG
	int RevLen = 0;
#endif
	i = PosQ - 1;
	j = PosT - 1;
	while (i >= 0 && j >= 0)
		{
		byte q = Q[i];
		byte t = T[j];
		RevScore += IntScoreMx_Mu[q][t];
		if (RevScore > BestRevScore)
			{
			BestRevScore = RevScore;
#if DEBUG
			Loi = i;
			Loj = j;
			RevLen = PosQ - i;
#endif
			}
		else if (RevScore + X < BestRevScore)
			break;
		--i;
		--j;
		}
	int BestScore = BestFwdScore + BestRevScore;
	EndTimer(MuXDrop);
#if DEBUG
	{
	Len = FwdLen + RevLen;
	int CheckScore = 0;
	int i = Loi;
	int j = Loj;
	for (int k = 0; k < Len; ++k)
		{
		asserta(i >= 0 && j >= 0 && i < LQ && j < LT);
		byte q = Q[i++];
		byte t = T[j++];
		CheckScore += IntScoreMx_Mu[q][t];
		}
	asserta(CheckScore == BestScore);
	}
#endif
	return BestScore;
	}

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
	int NewLo = NewScoreVec[MU_FILTER_KEEPN-1];
	s_QueryIdxToTargetIdxVec[QueryIdx] = NewTargetIdxVec;
	s_QueryIdxToScoreVec[QueryIdx] = NewScoreVec;
	s_QueryIdxToLoScore[QueryIdx] = NewLo;
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
	if (N == 0)
		{
		asserta(ScoreVec.empty());
		asserta(TargetIdxVec.empty());
		return;
		}
	int LoScore = s_QueryIdxToLoScore[QueryIdx];
	asserta(SIZE(TargetIdxVec) == N);
	map<uint, int> TargetIdxToScore;
	for (uint i = 0; i < N; ++i)
		{
		int Score = ScoreVec[i];
		asserta(Score >= LoScore);
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
	uint M = FullSize - 1;
	if (M > MU_FILTER_KEEPN-1)
		M = MU_FILTER_KEEPN-1;
	uint Lok = Order[M];
	int FullLoScore = FullScoreVec[Lok];

	for (uint k = 0; k < FullSize; ++k)
		{
		uint i = Order[k];
		assert(i < FullSize);
		uint TargetIdx = FullTargetIdxVec[i];
		int Score = FullScoreVec[i];
		map<uint, int>::const_iterator iter = TargetIdxToScore.find(TargetIdx);
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

static void CheckAllScoreVecs()
	{
	const uint QueryCount = SIZE(s_QueryIdxToScoreVec);
	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		CheckScoreVecs(QueryIdx);
	}
#endif // CHECK_SCORE_VECS

static void AddScore(uint QueryIdx, uint TargetIdx, int Score)
	{
	//brk(QueryIdx==0 && TargetIdx==23 && Score==66);
	s_DataLock.lock();
	++s_PveCount;
#if CHECK_SCORE_VECS
	s_QueryIdxToFullTargetIdxVec[QueryIdx].push_back(TargetIdx);
	s_QueryIdxToFullScoreVec[QueryIdx].push_back(Score);
#endif
	vector<int> &ScoreVec = s_QueryIdxToScoreVec[QueryIdx];
	int LoScore = s_QueryIdxToLoScore[QueryIdx];
	if (Score >= LoScore)
		{
		ScoreVec.push_back(Score);
		s_QueryIdxToTargetIdxVec[QueryIdx].push_back(TargetIdx);
		if (SIZE(ScoreVec) >= 2*MU_FILTER_KEEPN)
			TruncateVecs(QueryIdx);
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
					   const MuDex *ptrMD)
	{
	const DSSParams &Params = *ptrParams;
	FASTASeqSource &FSS = *ptrFSS;
	const MuDex &MD = *ptrMD;
	const string &PatternStr = Params.m_PatternStr;
	const uint QueryCount = s_ptrQueryDB->GetSeqCount();

	ObjMgr OM;

	SeqInfo *SI = OM.GetSeqInfo();
	vector<byte> MuLettersT;
	vector<uint> MuKmersT;
	int *QueryIdxToBestHSPScore = myalloc(int, QueryCount);
	for (;;)
		{
		bool Ok = FSS.GetNext(SI);
		if (!Ok)
			return;
		const uint LT = SI->m_L;
		const uint TargetIdx = SI->m_Index;
		const byte *ByteSeqT = SI->m_Seq;
		MuLettersT.clear();
		MuLettersT.reserve(LT);
		for (uint i = 0; i < LT; ++i)
			{
			byte c = ByteSeqT[i];
			byte Letter = g_CharToLetterMu[c];
			asserta(Letter < 36);
			MuLettersT.push_back(Letter);
			}
		const byte *ByteLettersT = MuLettersT.data();
		GetMuKmers(PatternStr, MuLettersT, MuKmersT);
		bool DoProgress = false;
		s_ProgressLock.lock();
		++s_TargetCount;
		if (1 || s_TargetCount%100 == 0)
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
			{
			double Pct = FSS.GetPctDone();
			Progress("%s chains scanned (%.1f%%)   \r",
					 IntToStr(s_TargetCount), Pct);
			}

		const uint KmerCountT = SIZE(MuKmersT);
		asserta(KmerCountT + 2 == LT);
		zero_array(QueryIdxToBestHSPScore, QueryCount);
		for (uint PosT = 0; PosT < KmerCountT; ++PosT)
			{
			uint KmerT = MuKmersT[PosT];
			uint RowSize = MD.GetRowSize(KmerT);
			const uint DataOffset = MD.GetDataOffset(KmerT);
			uint32_t QueryIdx;
			uint16_t PosQ;
			for (uint i = 0; i < RowSize; ++i)
				{
				MD.Get(DataOffset+i, QueryIdx, PosQ);
				assert(QueryIdx < QueryCount);
				assert(PosQ < s_ptrQueryDB->GetSeqLength(QueryIdx));
				uint LQ = s_QueryIdxToLength[QueryIdx];
				const byte *ByteLettersQ = s_QueryIdxToMuLetters[QueryIdx].data();
				int HSPScore = MuXDrop(
					ByteLettersQ, PosQ, LQ, 
					ByteLettersT, PosT, LT,
					X);
				if (HSPScore > MIN_HSP_SCORE)
					{
					if (HSPScore > QueryIdxToBestHSPScore[QueryIdx])
						QueryIdxToBestHSPScore[QueryIdx] = HSPScore;
					}
				}
			}
		for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
			{
			int Score = QueryIdxToBestHSPScore[QueryIdx];
			if (Score >= 30)
				{
#if CHECK_SCORE_VECSX
				asserta(GetRequestedThreadCount() == 1);
				CheckAllScoreVecs();
#endif
				AddScore(QueryIdx, TargetIdx, Score);
#if CHECK_SCORE_VECSX
				CheckAllScoreVecs();
#endif
				}
			}
		s_PairCount += QueryCount;
		}
	}

void cmd_mufilter()
	{
	asserta(optset_db);
	const string &QueryFN = g_Arg1;
	const string &DBFN = string(opt_db);
	FASTASeqSource FSS;
	FSS.Open(DBFN);

	if (!optset_output)
		Die("-output option required");
	FILE *fOut = CreateStdioFile(opt_output);

	DSSParams Params;
	Params.SetFromCmdLine(10000);
	const string &PatternStr = Params.m_PatternStr;
	asserta(PatternStr == "111");

	SeqDB QueryDB;
	QueryDB.FromFasta(QueryFN);
	const uint QueryCount = QueryDB.GetSeqCount();
	s_ptrQueryDB = &QueryDB;

	MuDex::Set_k(3);
	MuDex MD;
	MD.FromSeqDB(QueryDB);

	s_QueryIdxToScoreVec.resize(QueryCount);
	s_QueryIdxToTargetIdxVec.resize(QueryCount);
	s_QueryIdxToLoScore.resize(QueryCount, 0);
#if CHECK_SCORE_VECS
	s_QueryIdxToFullTargetIdxVec.resize(QueryCount);
	s_QueryIdxToFullScoreVec.resize(QueryCount);
#endif

	s_QueryIdxToMuLetters.resize(QueryCount);
	s_QueryIdxToMuKmers.resize(QueryCount);

	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		{
		ProgressStep(QueryIdx, QueryCount, "Indexing query");
		MuKmerFilter *MKF = new MuKmerFilter;
		MKF->SetParams(Params);
		const uint QL = QueryDB.GetSeqLength(QueryIdx);
		s_QueryIdxToLength.push_back(QL);
		const string &Q = QueryDB.GetSeq(QueryIdx);
		vector<byte> &MuLetters = s_QueryIdxToMuLetters[QueryIdx];
		vector<uint> &MuKmers = s_QueryIdxToMuKmers[QueryIdx];
		MuLetters.reserve(QL);
		MuKmers.reserve(QL);
		for (uint i = 0; i < QL; ++i)
			{
			byte c = Q[i];
			byte Letter = g_CharToLetterMu[c];
			asserta(Letter < 36);
			MuLetters.push_back(Letter);
			}
		GetMuKmers(PatternStr, MuLetters, MuKmers);
		}

	const uint ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody, ThreadIndex, &Params, &FSS, &MD);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

		Progress("%s chains scanned (100%%)    \n", IntToStr(s_TargetCount));

	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		TruncateVecs(QueryIdx);

	map<uint, vector<uint> > TargetIdxToQueryIdxs;
	vector<uint> TargetIdxs;
	for (uint QueryIdx = 0; QueryIdx < QueryCount; ++QueryIdx)
		{
		const vector<int> &ScoreVec = s_QueryIdxToScoreVec[QueryIdx];
		const vector<uint> &TargetIdxVec = s_QueryIdxToTargetIdxVec[QueryIdx];
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
	fprintf(fOut, "mufilter\t%u\n", TargetCount);
	for (uint k = 0; k < TargetCount; ++k)
		{
		uint TargetIdx = TargetIdxs[k];
		map<uint, vector<uint> >::const_iterator iter = TargetIdxToQueryIdxs.find(TargetIdx);
		asserta(iter != TargetIdxToQueryIdxs.end());
		const vector<uint> &QIdxs = iter->second;
		const uint K = SIZE(QIdxs);
		fprintf(fOut, "%u\t%u", TargetIdx, K);
		for (uint i = 0; i < K; ++i)
			fprintf(fOut, "\t%u", QIdxs[i]);
		fprintf(fOut, "\n");
		}
	CloseStdioFile(fOut);

#if CHECK_SCORE_VECS
	CheckAllScoreVecs();
	ProgressLog("CheckAllScoreVecs() ok\n");
#endif
	}
