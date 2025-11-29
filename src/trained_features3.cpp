#include "myutils.h"
#include "features.h"
#include "dss.h"

#pragma warning(disable:4305) // double -> float

float **g_ScoreMxs2[FEATURE_COUNT];
uint g_AlphaSizes2[FEATURE_COUNT];

void DSSParams::OverwriteUnweightedScoreMx(FEATURE F,
	vector<vector<float> > &ScoreMx)
	{
	asserta(g_ScoreMxs2[F] != 0);
	const uint AS = g_AlphaSizes2[F];
	asserta(SIZE(ScoreMx) == AS);
	for (uint i = 0; i < AS; ++i)
		{
		const vector<float> &Row = ScoreMx[i];
		asserta(SIZE(Row) == AS);
		for (uint j = 0; j < AS; ++j)
			g_ScoreMxs2[F][i][j] = ScoreMx[i][j];
		}
	}

static void SetFeatureScoreMx(FEATURE F, const float *mx, uint AS)
	{
	asserta(uint(F) < FEATURE_COUNT);
	asserta(AS == g_AlphaSizes2[F]);
	g_ScoreMxs2[F] = myalloc(float *, AS);
	for (uint i = 0; i < AS; ++i)
		{
		g_ScoreMxs2[F][i] = myalloc(float, AS);
		for (uint j = 0; j < AS; ++j)
			g_ScoreMxs2[F][i][j] = mx[AS*i + j];
		}
	}

static bool Init()
	{
#include "alphadata.inc"
	return true;
	}
static bool s_InitDone = Init();
