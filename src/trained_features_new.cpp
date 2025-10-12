#include "myutils.h"
#include "features.h"
#include "dss.h"

static void FreeFeature(FEATURE F);
static void AllocFeature(FEATURE F, uint AS);

// Do not make s_ScoreMxs public!
// DSSParams applies weights
static float **s_ScoreMxs[FEATURE_COUNT];
static float **s_FreqMxs[FEATURE_COUNT];
static float *s_FreqVecs[FEATURE_COUNT];
static uint s_AlphaSizes[FEATURE_COUNT];
static uint s_UndefLetters[FEATURE_COUNT];
static vector<vector<float> > s_BinTs;

static bool Init()
	{
	s_AlphaSizes[FEATURE_AA] = 20;
	s_AlphaSizes[FEATURE_Conf] = 16;
	s_AlphaSizes[FEATURE_NENConf] = 16;
	s_AlphaSizes[FEATURE_RENConf] = 16;
	s_AlphaSizes[FEATURE_PlusNENConf] = 16;
	s_AlphaSizes[FEATURE_MinusNENConf] = 16;
	s_AlphaSizes[FEATURE_Mu] = 36;

	s_BinTs.resize(FEATURE_COUNT);
	return true;
	};
static bool s_InitDone = Init();

static void AllocFeature(FEATURE F, uint AS)
	{
	s_AlphaSizes[F] = AS;
	s_FreqMxs[F] = myalloc(float *, AS);
	s_ScoreMxs[F] = myalloc(float *, AS);
	s_FreqVecs[F] = myalloc(float, AS);
	for (uint i = 0; i < AS; ++i)
		{
		s_FreqMxs[F][i] = myalloc(float, AS);
		s_ScoreMxs[F][i] = myalloc(float, AS);
		}
	}

static void FreeFeature(FEATURE F)
	{
	if (s_ScoreMxs[F] != 0)
		{
		uint AS = s_AlphaSizes[F];
		for (uint i = 0; i < AS; ++i)
			myfree(s_ScoreMxs[F][i]);
		myfree(s_ScoreMxs[F]);
		}
	if (s_FreqMxs[F] != 0)
		{
		uint AS = s_AlphaSizes[F];
		for (uint i = 0; i < AS; ++i)
			myfree(s_FreqMxs[F][i]);
		myfree(s_FreqMxs[F]);
		}
	if (s_FreqVecs[F] != 0)
		myfree(s_FreqVecs[F]);
	}

static void FreeMe()
	{
	for (uint F = 0; F < FEATURE_COUNT; ++F)
		FreeFeature(FEATURE(F));
	}

uint DSS::Get_PlusNENConf(uint Pos)
	{
	Die("Get_PlusNENConf");
	return UINT_MAX;
	}

uint DSS::Get_MinusNENConf(uint Pos)
	{
	Die("Get_MinusNENConf");
	return UINT_MAX;
	}

uint DSS::GetDefaultLetter(FEATURE F)
	{
	Die("GetDefaultLetter");
	assert(uint(F) < FEATURE_COUNT);
	return s_UndefLetters[uint(F)];
	}

const vector<float> &DSS::GetBinTs(FEATURE F)
	{
	Die("GetBinTs");
	assert(uint(F) < FEATURE_COUNT);
	return s_BinTs[uint(F)];
	}

void DSS::SetFeature(FEATURE F,
		const vector<float> &Freqs,
		const vector<vector<float> > &FreqMx,
		const vector<vector<float> > &ScoreMx,
		const vector<float> &BinTs,
		float UndefinedValue)
	{
	Die("SetFeature()");
#if 0
	asserta(uint(F) < FEATURE_COUNT);
	uint AS = SIZE(Freqs);
	asserta(SIZE(FreqMx) == AS);
	asserta(SIZE(ScoreMx) == AS);
	s_AlphaSizes[F] = AS;
	AllocFeature(F, AS);
	for (uint i = 0; i < AS; ++i)
		{
		s_FreqVecs[F][i] = Freqs[i];

		asserta(SIZE(FreqMx[i]) == AS);
		asserta(SIZE(ScoreMx[i]) == AS);

		for (uint j = 0; j < AS; ++j)
			{
			float Freq = FreqMx[i][j];
			float Score = ScoreMx[i][j];
			s_FreqMxs[F][i][j] = Freq;
			s_ScoreMxs[F][i][j] = Score;
			}
		}

	if (FeatureIsInt(F))
		asserta(BinTs.empty());
	else
		{
		uint n = GetBinThresholdCount(AS, UB);
		asserta(SIZE(BinTs) == n);
		s_BinTs[F].resize(n);
		if (!BinTs.empty())
			{
			for (uint i = 0; i < SIZE(BinTs); ++i)
				s_BinTs[F][i] = BinTs[i];
			}
		}
#endif
	}

const float *DSS::GetFreqVec(FEATURE F)
	{
	Die("GetFreqVec()");
	return s_FreqVecs[F];
	}

const float * const *DSS::GetFreqMx(FEATURE F)
	{
	Die("GetFreqMx()");
	return s_FreqMxs[F];
	}

const float * const *DSS::GetScoreMx(FEATURE F)
	{
	Die("GetScoreMx()");
	return s_ScoreMxs[F];
	}
