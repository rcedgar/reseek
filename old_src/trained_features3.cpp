#if 0
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

/***
* For experiment attempting to understand why new Conf subst
* matrices have slightly lower performance.
* up/down weights diagonal of Conf subst matrix, e.g. w=0.5 .. 1.5
***/
void SetConfDiagw(float w)
	{
	static const float S_ij[16*16] = {
	//     0        1        2        3        4        5        6        7        8        9       10       11       12       13       14       15
		1.09,   -4.33,   -4.42,   -4.62,   -3.48,   -4.51,   -4.33,   -1.02,   -2.24,   -2.89,    -4.7,  -0.296,   -2.15,   -4.22,   -4.23,   -3.98, // 0
	   -4.33,     1.7,    1.05,   0.374,    -1.5,   0.541,   -1.87,   -5.28,   -4.19,   -3.11,   0.617,   -3.62,   -5.47,   -1.11,   -2.52,    -4.4, // 1
	   -4.42,    1.05,    2.27,   -1.11,   -2.62, -0.0236,   -3.37,   -7.15,   -5.87,    -4.7,  -0.974,   -5.44,   -6.37,   -2.85,   -4.39,    -6.2, // 2
	   -4.62,   0.374,   -1.11,    2.11,   0.511,   0.634,   -1.11,   -2.31,   -1.77,   -1.72,  -0.101,   -1.16,   -2.23,   0.722,   0.594,   -2.47, // 3
	   -3.48,    -1.5,   -2.62,   0.511,    2.33,   0.638,   -1.35,   -3.42,    1.07,   -2.76,   -1.21,  -0.211,  -0.831,  -0.509,   0.252,   -2.19, // 4
	   -4.51,   0.541, -0.0236,   0.634,   0.638,    2.41,   -1.68,   -4.39,   -1.91,   -3.24,  -0.456,   -1.95,   -3.45,   -1.01,    -1.1,    -3.6, // 5
	   -4.33,   -1.87,   -3.37,   -1.11,   -1.35,   -1.68,    2.61,   -2.09,   0.341,   0.403,    1.31,   -3.01,    -3.1,  0.0272,   -0.52,   -1.05, // 6
	   -1.02,   -5.28,   -7.15,   -2.31,   -3.42,   -4.39,   -2.09,    2.84,   -2.17,  -0.772,   -2.85,  -0.978,    -1.7,  -0.809,  -0.163,  -0.343, // 7
	   -2.24,   -4.19,   -5.87,   -1.77,    1.07,   -1.91,   0.341,   -2.17,    2.47,  -0.426,   -1.35,  -0.277,   0.735,    -1.3,   0.151, -0.0869, // 8
	   -2.89,   -3.11,    -4.7,   -1.72,   -2.76,   -3.24,   0.403,  -0.772,  -0.426,    2.76,   -1.24,    -2.5,   -2.96,    1.07,   -1.06,    1.52, // 9
		-4.7,   0.617,  -0.974,  -0.101,   -1.21,  -0.456,    1.31,   -2.85,   -1.35,   -1.24,     2.5,   -3.08,   -4.28,       0,   -1.25,    -2.3, // 10
	  -0.296,   -3.62,   -5.44,   -1.16,  -0.211,   -1.95,   -3.01,  -0.978,  -0.277,    -2.5,   -3.08,     2.3,   0.768,  -0.562,  0.0675,   -2.14, // 11
	   -2.15,   -5.47,   -6.37,   -2.23,  -0.831,   -3.45,    -3.1,    -1.7,   0.735,   -2.96,   -4.28,   0.768,    2.98,   -1.62,   0.935,   0.368, // 12
	   -4.22,   -1.11,   -2.85,   0.722,  -0.509,   -1.01,  0.0272,  -0.809,    -1.3,    1.07,       0,  -0.562,   -1.62,     2.6,   0.254,   0.191, // 13
	   -4.23,   -2.52,   -4.39,   0.594,   0.252,    -1.1,   -0.52,  -0.163,   0.151,   -1.06,   -1.25,  0.0675,   0.935,   0.254,    2.63, -0.0195, // 14
	   -3.98,    -4.4,    -6.2,   -2.47,   -2.19,    -3.6,   -1.05,  -0.343, -0.0869,    1.52,    -2.3,   -2.14,   0.368,   0.191, -0.0195,     3.1, // 15
	};
	for (uint i = 0; i < 16; ++i)
		g_ScoreMxs2[FEATURE_Conf][i][i] = S_ij[16*i + i]*w;
	}
#endif // 0