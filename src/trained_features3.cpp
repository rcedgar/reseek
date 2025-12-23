#include "myutils.h"
#include "features.h"
#include "dss.h"

#pragma warning(disable:4305) // double -> float
#pragma warning(disable:4244) // int -> float

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

void DSSParams::CreateFeatureScoreMx(FEATURE F, 
	const vector<vector<float> > &ScoreMx)
	{
	const uint AS = SIZE(ScoreMx);
	asserta(uint(F) < FEATURE_COUNT);
	g_AlphaSizes2[F] = AS;
	g_ScoreMxs2[F] = myalloc(float *, AS);
	for (uint i = 0; i < AS; ++i)
		{
		g_ScoreMxs2[F][i] = myalloc(float, AS);
		for (uint j = 0; j < AS; ++j)
			g_ScoreMxs2[F][i][j] = ScoreMx[i][j];
		}
	}

static void SetFeatureScoreMx(FEATURE F, const float *mx, uint AS)
	{
	asserta(uint(F) < FEATURE_COUNT);
	g_AlphaSizes2[F] = AS;
	g_ScoreMxs2[F] = myalloc(float *, AS);
	for (uint i = 0; i < AS; ++i)
		{
		g_ScoreMxs2[F][i] = myalloc(float, AS);
		for (uint j = 0; j < AS; ++j)
			g_ScoreMxs2[F][i][j] = mx[AS*i + j];
		}
	}

static void SetPlaceholderScoreMx(FEATURE F, uint AS)
	{
	asserta(uint(F) < FEATURE_COUNT);
	g_AlphaSizes2[F] = AS;
	g_ScoreMxs2[F] = myalloc(float *, AS);
	for (uint i = 0; i < AS; ++i)
		{
		g_ScoreMxs2[F][i] = myalloc(float, AS);
		for (uint j = 0; j < AS; ++j)
			g_ScoreMxs2[F][i][j] = (i == j ? 1 : 0);
		}
	}

static bool Init()
	{
	// SetPlaceholderScoreMx(FEATURE_MENDist4b, 4);

#include "alphadata.h"
	return true;
	}
static bool s_InitDone = Init();

void SetBLAST_B62()
	{
	extern int Blosum62_int[20][20];
	for (uint i = 0; i < 20; ++i)
		for (uint j = 0; j < 20; ++j)
			g_ScoreMxs2[FEATURE_B62][i][j] = Blosum62_int[i][j];
	}
