#include "myutils.h"
#include "chainer.h"
#include "sort.h"

#define TRACE	0
#define TEST	0

static const float MINUS_INFINITY = -9e9f;

// Ties: Los before His
static int CmpBPs(const void *vpBP1, const void *vpBP2)
	{
	const BPData *BP1 = (const BPData *) vpBP1;
	const BPData *BP2 = (const BPData *) vpBP2;

	if (BP1->Pos < BP2->Pos)
		return -1;
	else if (BP1->Pos > BP2->Pos)
		return 1;
	assert(BP1->Pos == BP2->Pos);
	if (BP1->IsLo != BP2->IsLo)
		{
		if (BP1->IsLo && !BP2->IsLo)
			return -1;
		else
			return 1;
		}
	return 0;
	}

void Chainer::Clear()
	{
	if (m_BPs != 0) { myfree(m_BPs); m_BPs = 0; }
	if (m_ChainScores != 0) { myfree(m_ChainScores); m_ChainScores = 0; }
	}

float Chainer::Chain(const vector<uint> &Los, const vector<uint> &His,
	  const vector<float> &Scores, vector<uint> &Idxs)
	{
	Idxs.clear();
	const uint N = SIZE(Los);
	assert(SIZE(His) == N);
	assert(SIZE(Los) == N);
	if (N == 0)
		return 0;

#if	TRACE
	Log("Chainer::Chain(N=%u)\n", N);
#endif
	Clear();
	m_BPs = myalloc(BPData, 2*N);

	for (uint i = 0; i < N; ++i)
		{
		uint Lo = Los[i];
		uint Hi = His[i];
		asserta(Hi >= Lo);

		m_BPs[2*i].Index = i;
		m_BPs[2*i].IsLo = true;
		m_BPs[2*i].Pos = Lo;

		m_BPs[2*i+1].Index = i;
		m_BPs[2*i+1].IsLo = false;
		m_BPs[2*i+1].Pos = Hi;
		}

#if	0 // TRACE
	{
	Log("BPs:\n");
	Log("    Pos    Index  LH   Score\n");
	Log("-------  -------  --  ------\n");
	for (uint i = 0; i < 2*N; ++i)
		{
		const BPData &BP = m_BPs[i];
		Log("%7u", BP.Pos);
		Log("  %7u", BP.Index);
		Log("  %s", BP.IsLo ? "Lo" : "Hi");
		Log("  %6.1f", Scores[BP.Index]);
		Log("\n");
		}
	}
#endif

	//SortBPVecInPlace(m_BPs, 2*N);
	qsort(m_BPs, 2*N, sizeof(BPData), CmpBPs);

#if	TRACE
	{
	Log("Sorted BPs:\n");
	Log("    Pos    Index  LH   Score\n");
	Log("-------  -------  --  ------\n");
	for (uint i = 0; i < 2*N; ++i)
		{
		const BPData &BP = m_BPs[i];
		Log("%7u", BP.Pos);
		Log("  %7u", BP.Index);
		Log("  %s", BP.IsLo ? "Lo" : "Hi");
		Log("  %6.1f", Scores[BP.Index]);
		Log("\n");
		}
	}
#endif

	m_TB = myalloc(uint, N);
	m_ChainScores = myalloc(float, N);

	asserta(m_BPs[0].IsLo);
	uint Index0 = m_BPs[0].Index;
	uint BestChainEnd = UINT_MAX;
	m_TB[0] = UINT_MAX;

	m_ChainScores[0] = MINUS_INFINITY;

	for (uint i = 0; i < 2*N; ++i)
		{
		const BPData &BP = m_BPs[i];

		assert(BP.Index < N);
		float Score = Scores[BP.Index];

		if (BP.IsLo)
			{
			m_TB[BP.Index] = BestChainEnd;
			if (BestChainEnd == UINT_MAX)
				m_ChainScores[BP.Index] = Score;
			else
				m_ChainScores[BP.Index] = m_ChainScores[BestChainEnd] + Score;
			}
		else
			{
			if (BestChainEnd == UINT_MAX || m_ChainScores[BP.Index] > m_ChainScores[BestChainEnd])
				BestChainEnd = BP.Index;
			}
		}

	asserta(BestChainEnd < N);

#if	TRACE
	{
	Log("\n");
	Log("BestChainEnd %u, Score %.1f\n", BestChainEnd, m_ChainScores[BestChainEnd]);
	Log("Index  ChainScore     TB\n");
	Log("-----  ----------  -----\n");
	for (uint i = 0; i < N; ++i)
		{
		Log("%5u", i);
		float Score = m_ChainScores[i];
		if (Score == MINUS_INFINITY)
			Log("  %10.10s", "*");
		else
			Log("  %10.1f", Score);
		uint t = m_TB[i];
		if (t == UINT_MAX)
			Log("  %5.5s", "*");
		else
			Log("  %5u", t);
		Log("\n");
		}
	}
#endif

	uint Index = BestChainEnd;
	float Score = 0;
	for (;;)
		{
		assert(SIZE(Idxs) < N);
		Score += Scores[Index];
		Idxs.push_back(Index);
		assert(Index < N);
		Index = m_TB[Index];
		if (Index == UINT_MAX)
			break;
		}

#if	TRACE
	{
	Log("\n");
	Log("Chain:\n");
	Log("Index     Lo     Hi   Score\n");
	Log("-----  -----  -----  ------\n");
	float Sum = 0.0;
	const uint ChainLength = SIZE(Idxs);
	for (uint i = 0; i < ChainLength; ++i)
		{
		uint Index = Idxs[i];
		asserta(Index < N);
		Log("%5u", Index);
		Log("  %5u", Los[Index]);
		Log("  %5u", His[Index]);
		Log("  %6.1f", Scores[Index]);
		Sum += Scores[Index];
		Log("\n");
		}
	Log("Sum %.1f\n", Sum);
	}
#endif
	return Score;
	}

float Chainer::GetChainScore(const uint *Los, const uint *His,
  const float *Scores, uint N, const vector<uint> &Idxs)
	{
	const uint ChainLength = SIZE(Idxs);
	float Sum = 0.0;
	for (uint i = 0; i < ChainLength; ++i)
		{
		uint Index = Idxs[i];
		assert(Index < N);
		Sum += Scores[Index];
		}
	return Sum;
	}

#if	TEST
const uint MaxTries = 100;
const uint MinCount = 1;
const uint MaxCount = 8;
const uint MinLen = 1;
const uint MaxLen = 100;
const uint MaxPos = 100;
const uint MinScore = 1;
const uint MaxScore = 100;
const uint RandSeed = 0;

static void GetRandomLoHi(uint MaxPos, uint MinLen, uint MaxLen,
  uint MinScore, uint MaxScore, uint &Lo, uint &Hi, float &Score)
	{
	asserta(MinLen <= MaxLen);
	asserta(MinScore <= MaxScore);

	Lo = uint(rand()%MaxPos);
	uint Length = MinLen + uint(rand()%(MaxLen - MinLen + 1));
	Hi = Lo + Length - 1;
	Score = float(MinScore + uint(rand()%(MaxScore - MinScore + 1)));
	}

static uint GetRandomLoHis(
  uint MinCount, uint MaxCount,
  uint MaxPos,
  uint MinLen, uint MaxLen,
  uint MinScore, uint MaxScore,
  uint *Los, uint *His, float *Scores)
	{
	asserta(MinCount <= MaxCount);
	uint Count = MinCount + uint(rand()%(MaxCount - MinCount + 1));
	for (uint i = 0; i < Count; ++i)
		GetRandomLoHi(MaxPos, MinLen, MaxLen, MinScore, MaxScore,
		  Los[i], His[i], Scores[i]);
	return Count;
	}

void cmd_test()
	{
	srand(RandSeed);

	Chainer C;

	uint *Los = myalloc(uint, MaxCount);
	uint *His = myalloc(uint, MaxCount);
	float *Scores = myalloc(float, MaxCount);

	for (uint Try = 0; Try < MaxTries; ++Try)
		{
		ProgressStep(Try, MaxTries, "Testing");
		uint N = GetRandomLoHis(MinCount, MaxCount, MaxPos, MinLen, MaxLen, MinScore, MaxScore,
		  Los, His, Scores);

		vector<uint> Idxs;
		C.Chain(Los, His, Scores, N, Idxs);
		float Score = Chainer::GetChainScore(Los, His, Scores, N, Idxs);
		const uint ChainLength = SIZE(Idxs);
		Log("N %u, chain %u, Score %.1f\n", N, ChainLength, Score);
		}
	}
#endif // TEST
