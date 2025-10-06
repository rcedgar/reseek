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
static uint s_DefaultLetters[FEATURE_COUNT];
static UNDEF_BINNING s_UBs[FEATURE_COUNT];
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

//uint DSS::GetAlphaSize(FEATURE F)
//	{
//	assert(uint(F) < FEATURE_COUNT);
//	return s_AlphaSizes[uint(F)];
//	}

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
	assert(uint(F) < FEATURE_COUNT);
	return s_DefaultLetters[uint(F)];
	}

UNDEF_BINNING DSS::GetUB(FEATURE F)
	{
	assert(uint(F) < FEATURE_COUNT);
	return s_UBs[uint(F)];
	}

const vector<float> &DSS::GetBinTs(FEATURE F)
	{
	assert(uint(F) < FEATURE_COUNT);
	return s_BinTs[uint(F)];
	}

void DSS::SetFeature(FEATURE F, UNDEF_BINNING UB,
		const vector<float> &Freqs,
		const vector<vector<float> > &FreqMx,
		const vector<vector<float> > &ScoreMx,
		const vector<float> &BinTs)
	{
	asserta(uint(F) < FEATURE_COUNT);
	uint AS = SIZE(Freqs);
	asserta(SIZE(FreqMx) == AS);
	asserta(SIZE(ScoreMx) == AS);
	s_AlphaSizes[F] = AS;
	s_UBs[F] = UB;
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
	}

const float *DSS::GetFreqVec(FEATURE F)
	{
	return s_FreqVecs[F];
	}

const float * const *DSS::GetFreqMx(FEATURE F)
	{
	return s_FreqMxs[F];
	}

const float * const *DSS::GetScoreMx(FEATURE F)
	{
	return s_ScoreMxs[F];
	}

//uint DSS::ValueToInt_Feature(FEATURE F, float Value)
//	{
//	assert(uint(F) < FEATURE_COUNT);
//	uint AS = s_AlphaSizes[F];
//	assert(AS > 0);
//	UNDEF_BINNING UB = s_UBs[F];
//	const vector<float> &BinTs = s_BinTs[F];
//	uint DefaultLetter = s_DefaultLetters[F];
//	uint Letter = ValueToInt(Value, UB, AS, BinTs, DefaultLetter);
//	return Letter;
//	}
