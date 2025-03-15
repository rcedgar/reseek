#include "myutils.h"
#include "features.h"
#include "dss.h"

static float **g_ScoreMxs3[FEATURE_COUNT];
static float **g_FreqMxs3[FEATURE_COUNT];
static float *g_FreqVecs3[FEATURE_COUNT];
static float *g_BinTs3[FEATURE_COUNT];
static uint g_AlphaSizes3[FEATURE_COUNT];

const float *DSS::GetNewFreqVec(FEATURE F)
	{
	return g_FreqVecs3[F];
	}

const float * const *DSS::GetNewFreqMx(FEATURE F)
	{
	return g_FreqMxs3[F];
	}

const float * const *DSS::GetNewScoreMx(FEATURE F)
	{
	return g_ScoreMxs3[F];
	}

static void AllocNewFeature(FEATURE F, uint AS)
	{
	g_AlphaSizes3[F] = AS;
	g_FreqMxs3[F] = myalloc(float *, AS);
	g_ScoreMxs3[F] = myalloc(float *, AS);
	g_FreqVecs3[F] = myalloc(float, AS);
	g_BinTs3[F] = myalloc(float, AS);
	for (uint i = 0; i < AS; ++i)
		{
		g_FreqMxs3[F][i] = myalloc(float, AS);
		g_ScoreMxs3[F][i] = myalloc(float, AS);
		}
	}

static void FreeNewFeature(FEATURE F)
	{
	if (g_ScoreMxs3[F] != 0)
		{
		uint AS = g_AlphaSizes3[F];
		for (uint i = 0; i < AS; ++i)
			myfree(g_ScoreMxs3[F][i]);
		myfree(g_ScoreMxs3[F]);
		}
	if (g_FreqMxs3[F] != 0)
		{
		uint AS = g_AlphaSizes3[F];
		for (uint i = 0; i < AS; ++i)
			myfree(g_FreqMxs3[F][i]);
		myfree(g_FreqMxs3[F]);
		}
	if (g_FreqVecs3[F] != 0)
		myfree(g_FreqVecs3[F]);
	if (g_BinTs3[F] != 0)
		myfree(g_BinTs3[F]);
	}

void DSS::SetNewTrainFeature(FEATURE F,
		const vector<float> &Freqs,
		const vector<vector<float> > &FreqMx,
		const vector<vector<float> > &ScoreMx,
		const vector<float> &BinTs)
	{
	asserta(uint(F) < FEATURE_COUNT);
	uint AS = SIZE(Freqs);
	asserta(SIZE(FreqMx) == AS);
	asserta(SIZE(ScoreMx) == AS);
	g_AlphaSizes3[uint(F)] = AS;
	AllocNewFeature(F, AS);
	for (uint i = 0; i < AS; ++i)
		{
		g_FreqVecs3[F][i] = Freqs[i];
		g_BinTs3[F][i] = BinTs[i];

		asserta(SIZE(FreqMx[i]) == AS);
		asserta(SIZE(ScoreMx[i]) == AS);

		for (uint j = 0; j < AS; ++j)
			{
			(*g_FreqMxs3)[i][j] = FreqMx[i][j];
			(*g_ScoreMxs3)[i][j] = ScoreMx[i][j];
			}
		}
	}
