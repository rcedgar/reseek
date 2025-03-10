#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"

/***
  -sscluster d:/a/res/dave_grant/scop40/scop40.tm0.6_0.8.fa2
  -train_cal d:/int/scop40/out/domains_scop.cal 
  -myss3 Y 
  -randseed 6 
  -k 16 
  -n 100000 
  -log k16.seed6.ss3Y 
***/

const uint M = 9;
const uint K = 16;
static vector<int> ivalues;
static vector<int> jvalues;
static vector<vector<double> > Means;

static void InitLetter(uint Letter, uint Count,
  double x0, double x1, double x2, double x3, double x4,
  double x5, double x6, double x7, double x8)
	{
	asserta(Letter < K);
	Means[Letter][0] = x0;
	Means[Letter][1] = x1;
	Means[Letter][2] = x2;
	Means[Letter][3] = x3;
	Means[Letter][4] = x4;
	Means[Letter][5] = x5;
	Means[Letter][6] = x6;
	Means[Letter][7] = x7;
	Means[Letter][8] = x8;
	}

static bool Init()
	{
	Means.resize(K);
	for (uint k = 0; k < K; ++k)
		Means[k].resize(M, DBL_MAX);

	for (int i = -2; i <= 2; ++i)
		for (int j = i+1; j <= 2; ++j)
			{
			int minij = min(i, j);
			int maxij = max(i, j);
			if (maxij - minij == 1)
				continue;
			ivalues.push_back(minij);
			jvalues.push_back(maxij);
			}

	ivalues.push_back(-3);
	jvalues.push_back(3);

	ivalues.push_back(0);
	jvalues.push_back(3);

	ivalues.push_back(-3);
	jvalues.push_back(0);

	asserta(SIZE(ivalues) == M);
	asserta(SIZE(jvalues) == M);

#define SSKMEAN(Letter, Count, x0, x1, x2, x3, x4, x5, x6, x7, x8) \
	InitLetter(Letter, Count, x0, x1, x2, x3, x4, x5, x6, x7, x8);

//                             -2,0        -2,1        -2,2        -1,1        -1,2         0,2        -3,3         0,3        -3,0
SSKMEAN(  0,      25988,      5.466,      5.218,      6.295,      5.472,      5.231,      5.478,      9.957,      5.326,      5.203);
SSKMEAN(  1,      10773,      6.561,      9.619,       12.5,      6.584,       9.54,      6.503,      17.84,      9.217,      9.588);
SSKMEAN(  2,       8482,      6.872,      10.26,      13.44,      6.863,      10.24,      6.839,      19.66,      10.12,      10.16);
SSKMEAN(  3,       5319,      6.003,      8.082,      10.41,      6.348,      8.912,      6.544,      15.41,      9.469,      8.833);
SSKMEAN(  4,       5188,      5.795,      8.276,      10.66,      6.402,      9.297,      6.554,      12.19,      9.211,      5.914);
SSKMEAN(  5,       4953,      6.581,      9.627,      12.72,      6.623,      9.836,      6.776,      15.52,      9.617,      7.687);
SSKMEAN(  6,       4730,      6.506,      9.369,      11.06,      6.583,      8.264,      5.669,      11.32,      5.711,      8.937);
SSKMEAN(  7,       4459,      5.573,      5.537,      6.667,      5.473,      5.418,      5.498,      11.02,      5.621,      8.287);
SSKMEAN(  8,       4421,      5.679,      7.569,      9.335,      6.127,      8.475,      5.949,      8.031,      7.051,      5.713);
SSKMEAN(  9,       4279,      6.423,      8.457,      7.479,      5.695,      5.693,      5.486,      10.18,      5.738,      9.094);
SSKMEAN( 10,       4262,      6.636,      9.453,      11.96,      6.413,      9.031,      5.945,      14.73,      6.715,      9.583);
SSKMEAN( 11,       4196,      5.474,       5.71,      7.441,      5.583,      6.551,      6.021,         12,      8.539,      5.595);
SSKMEAN( 12,       3582,      5.574,      5.368,      6.191,      5.682,      7.545,      6.115,      6.775,      8.909,      6.277);
SSKMEAN( 13,       3366,      6.387,      8.708,      8.938,      5.809,      6.408,      6.041,      13.26,      8.584,      9.058);
SSKMEAN( 14,       3290,      5.867,      6.271,      8.535,      5.848,      8.319,      6.471,      10.82,      9.066,      8.534);
SSKMEAN( 15,       2712,      6.515,      7.566,      6.319,      5.675,      5.571,      5.653,      6.643,      7.279,      8.784);

#undef SKMEAN
	return true;
	}
static bool g_InitDone = Init();

static double Conf[16][16] = {
//                   A           C           D           E           F           G           H           I           K           L           M           N           P           Q           R           S
/*   A */ {      2.731,     0.5676,    -0.9645,     -1.992,    0.01539,     -2.716,     0.7432,     0.8708,     -1.487,     -1.101,     -0.928,     -1.114,    -0.8759,     0.2462,    -0.7739,    -0.9635}, // A
/*   C */ {     0.5676,      2.548,     0.3908,     -1.252,     0.4374,     -2.252,     0.2289,     0.4025,    -0.4552,    -0.1668,    -0.1715,      -1.91,    0.06908,     0.1485,     0.5235,     -0.377}, // C
/*   D */ {    -0.9645,     0.3908,      1.967,     0.2477,     0.4142,    -0.9734,     -1.098,    -0.7217,    0.04012,     -1.111,     -0.822,     -2.467,    -0.4636,    -0.6726,     0.6847,     0.4521}, // D
/*   E */ {     -1.992,     -1.252,     0.2477,      1.411,    -0.8166,     0.6886,     -1.296,     -1.857,     0.4187,     -2.204,     -1.359,     -3.229,      -0.93,     -1.693,    -0.6826,     0.3617}, // E
/*   F */ {    0.01539,     0.4374,     0.4142,    -0.8166,      2.041,     -1.431,    -0.6386,     0.6659,    -0.6156,     -1.142,    -0.9187,     -1.487,     -0.394,     0.0874,    0.01761,     0.5354}, // F
/*   G */ {     -2.716,     -2.252,    -0.9734,     0.6886,     -1.431,      1.847,     -2.477,     -2.959,    -0.7155,     -3.274,      -2.38,     -4.328,     -1.748,     -2.782,      -1.81,    -0.1411}, // G
/*   H */ {     0.7432,     0.2289,     -1.098,     -1.296,    -0.6386,     -2.477,      2.946,     0.4304,     -0.507,    -0.1757,      1.098,     -1.874,   -0.09558,    -0.9906,    -0.1277,     -1.373}, // H
/*   I */ {     0.8708,     0.4025,    -0.7217,     -1.857,     0.6659,     -2.959,     0.4304,      2.196,    -0.7989,    -0.5996,   0.002434,    -0.8155,     0.2228,     0.1478,    -0.4203,    -0.8076}, // I
/*   K */ {     -1.487,    -0.4552,    0.04012,     0.4187,    -0.6156,    -0.7155,     -0.507,    -0.7989,       2.16,     -1.203,     -0.428,     -2.178,     0.8158,     -1.184,     0.1879,    -0.3276}, // K
/*   L */ {     -1.101,    -0.1668,     -1.111,     -2.204,     -1.142,     -3.274,    -0.1757,    -0.5996,     -1.203,      2.533,    0.07438,    -0.2165,     -0.637,    -0.1676,    -0.3425,     -1.813}, // L
/*   M */ {     -0.928,    -0.1715,     -0.822,     -1.359,    -0.9187,      -2.38,      1.098,   0.002434,     -0.428,    0.07438,       2.45,     -1.016,      0.421,    -0.8731,     0.6288,     -1.435}, // M
/*   N */ {     -1.114,      -1.91,     -2.467,     -3.229,     -1.487,     -4.328,     -1.874,    -0.8155,     -2.178,    -0.2165,     -1.016,      1.116,     -1.557,     0.1198,      -1.91,     -2.535}, // N
/*   P */ {    -0.8759,    0.06908,    -0.4636,      -0.93,     -0.394,     -1.748,   -0.09558,     0.2228,     0.8158,     -0.637,      0.421,     -1.557,      2.271,    -0.8715,     0.2738,     -0.815}, // P
/*   Q */ {     0.2462,     0.1485,    -0.6726,     -1.693,     0.0874,     -2.782,    -0.9906,     0.1478,     -1.184,    -0.1676,    -0.8731,     0.1198,    -0.8715,      2.191,    -0.1029,    -0.8234}, // Q
/*   R */ {    -0.7739,     0.5235,     0.6847,    -0.6826,    0.01761,      -1.81,    -0.1277,    -0.4203,     0.1879,    -0.3425,     0.6288,      -1.91,     0.2738,    -0.1029,      2.458,    -0.4951}, // R
/*   S */ {    -0.9635,     -0.377,     0.4521,     0.3617,     0.5354,    -0.1411,     -1.373,    -0.8076,    -0.3276,     -1.813,     -1.435,     -2.535,     -0.815,    -0.8234,    -0.4951,      2.063}, // S
};

static double GetDist(
  const vector<double> &v1,
  const vector<double> &v2)
	{
	const uint n = SIZE(v1);
	asserta(SIZE(v2) == n);
	double Sum2 = 0;
	for (uint i = 0; i < n; ++i)
		{
		double diff = v1[i] - v2[i];
		Sum2 += diff*diff;
		}
	return sqrt(Sum2);
	}

static uint GetConfLetter(const vector<double> &v)
	{
	asserta(SIZE(v) == M);
	double MinDist = DBL_MAX;
	uint BestCluster = WILDCARD;
	for (uint k = 0; k < K; ++k)
		{
		double d = GetDist(v, Means[k]);
		if (k == 0 || d < MinDist)
			{
			BestCluster = k;
			MinDist = d;
			}
		}
	return BestCluster;
	}

static void Getv(const PDBChain &Chain, uint Pos,
  vector<double> &v)
	{
	v.clear();
	const uint L = Chain.GetSeqLength();
	if (Pos < 3 || Pos + 3 >= int(L))
		return;

	for (uint m = 0; m < M; ++m)
		{
		int i = ivalues[m];
		int j = jvalues[m];
		double d = Chain.GetDist(Pos+i, Pos+j);
		v.push_back(d);
		}
	asserta(SIZE(v) == M);
	}

uint DSS::Get_Conf(uint Pos)
	{
	vector<double> v;
	Getv(*m_Chain, Pos, v);
	if (v.empty())
		return WILDCARD;
	uint Letter = GetConfLetter(v);
	return Letter;
	}

uint DSS::Get_NENConf(uint Pos)
	{
	SetNENs();
	vector<double> v;
	uint NEN = GetNEN(Pos);
	if (NEN == UINT_MAX)
		return WILDCARD;

	Getv(*m_Chain, NEN, v);
	if (v.empty())
		return WILDCARD;

	uint Letter = GetConfLetter(v);
	if (Letter == UINT_MAX)
		return WILDCARD;
	return Letter;
	}

uint DSS::Get_PlusNENConf(uint Pos)
	{
	SetNENs();
	vector<double> v;
	uint NEN = GetPlusNEN(Pos);
	if (NEN == UINT_MAX)
		return WILDCARD;

	Getv(*m_Chain, NEN, v);
	if (v.empty())
		return WILDCARD;

	uint Letter = GetConfLetter(v);
	if (Letter == UINT_MAX)
		return WILDCARD;
	return Letter;
	}

uint DSS::Get_MinusNENConf(uint Pos)
	{
	SetNENs();
	vector<double> v;
	uint NEN = GetMinusNEN(Pos);
	if (NEN == UINT_MAX)
		return WILDCARD;

	Getv(*m_Chain, NEN, v);
	if (v.empty())
		return WILDCARD;

	uint Letter = GetConfLetter(v);
	if (Letter == UINT_MAX)
		return WILDCARD;
	return Letter;
	}

uint DSS::Get_RENConf(uint Pos)
	{
	SetNENs();
	vector<double> v;
	uint NEN = GetREN(Pos);
	if (NEN == UINT_MAX)
		return WILDCARD;

	Getv(*m_Chain, NEN, v);
	if (v.empty())
		return WILDCARD;

	uint Letter = GetConfLetter(v);
	if (Letter == UINT_MAX)
		return WILDCARD;
	asserta(Letter < 16);
	return Letter;
	}
