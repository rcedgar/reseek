#include "myutils.h"
#include "flat_chain.h"
#include "flat_distmx.h"
#include "flat_aligner.h"
#include "flat_helpers.h"
#include "sort.h"

static const float tx = 1.25;//TODO param
//static float s_DALI_Theta = 0;
float g_DALI_D = 20.0f;
float g_DALI_d0 = 0.2f;
float g_DALI_Theta = 1.0f;

static const int TBLSZ = 100;
static double *WeightLookup;

/***
DaliLite v5
comparemodules.f, line 1436
===========================
		enveloperadius=20.0
		x=1/(enveloperadius*enveloperadius)
		do i=0,100
				wght(i)=exp(-x*i*i)
		end do
***/
static double Weight(double y)
	{
	int iy = int(y+0.5);
	if (iy < 0)
		iy = 0;
	if (iy >= TBLSZ)
		iy = TBLSZ-1;
	double w2 = WeightLookup[iy];
	return w2;
	}

static double Weight_NoLookup(double y)
	{
	//const double D = 20.0;
	const double x = 1.0 / (g_DALI_D * g_DALI_D);
	double w = exp(-x * y * y);
	return w;
	}

static void FreeMe();

static bool InitWeightLookup()
	{
	atexit(FreeMe);
	WeightLookup = myalloc(double, TBLSZ);
	for (int i = 0; i < TBLSZ; ++i)
		{
		double y = double(i);
		double w = Weight_NoLookup(y);
		WeightLookup[i] = w;
		}
	return true;
	};
static bool InitWeightLookupDone = InitWeightLookup();

static void FreeMe()
	{
	myfree(WeightLookup);
	}

/***
comparemodules.f, line 1397
  a, b are integer distances in units of 1/10 Angstrom,
  so multiply by 10 to get Angstroms.
===========================
		function dpscorefun(a,b) result(s)
		implicit none
		include 'parsizes.for'
		real s
		integer*2 a,b
c
		real x,y,d0
		logical lela
		parameter(lela=.true.)
		parameter(d0=0.20)
c !!!   elastic uses weights !!!
		x=float(abs(a-b))/10
		if(lela) then
				y=float(a+b)/20
				if(y.gt.100) then
						s=0.0
				else
						if(y.gt.0) then
						  s=wght(nint(y))*(d0-x/y)
						else
						  s=wght(nint(y))*d0
						end if
				end if
		end if
***/
double DALI_dpscorefun(double a, double b)
	{
	double Score = 0;
	double diff = fabs(a - b);
	double mean = (a + b) / 2;
	double ratio = diff/mean;
	double w = Weight(mean);
	if (mean > 100)
		Score = 0;
	else
		{
		if (mean > 0)
			Score = w*(g_DALI_d0 - ratio);
		else
			Score = w*g_DALI_d0;
		}
	return Score;
	}

float flat_get_dali4(
	const uint *posQs, uint LQ, 
	const uint *posTs, uint LT, uint nmatch,
	const sid_t *distmxQ, const sid_t *distmxT)
	{
	const uint M = flat_params::m_distmx_bandwidth;

	extern float g_DALI_Theta;
	float score = g_DALI_Theta*nmatch;
	for (uint coli = 0; coli < nmatch; ++coli)
		{
		uint posQi = posQs[coli];
		uint posTi = posTs[coli];
		assert(posQi != UINT_MAX);
		assert(posTi != UINT_MAX);
		assert(posQi < LQ);
		assert(posTi < LT);

		for (uint colj = coli + 1; colj < nmatch; ++colj)
			{
			uint posQj = posQs[colj];
			uint posTj = posTs[colj];
			assert(posQj != UINT_MAX);
			assert(posTj != UINT_MAX);
			assert(posQj < LQ);
			assert(posTj < LT);

			int diffij_Q = abs(int(posQi) - int(posQj));
			int diffij_T = abs(int(posTi) - int(posTj));
			if (diffij_Q > int(M) || diffij_T > int(M))
				continue;
			uint kQ = banded_ij_to_k(posQi, posQj);
			uint kT = banded_ij_to_k(posTi, posTj);
			sid_t sid_dQ_squared = distmxQ[kQ];
			sid_t sid_dT_squared = distmxT[kT];
			float dQ_squared = sid2dist2(sid_dQ_squared);
			float dT_squared = sid2dist2(sid_dT_squared);
			float dQ = sqrtf(dQ_squared);
			float dT = sqrtf(dT_squared);

			double DALI_dpscorefun(double a, double b);
			score += (float) DALI_dpscorefun(dQ, dT);
			}
		}
	return score;
	}

float flat_get_dali2(
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const vector<uint> &posQs,
	const vector<uint> &posTs,
	const sid_t *distmxQ,
	const sid_t *distmxT)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint ncol = uint(posQs.size());

	extern float g_DALI_Theta;
	float score = g_DALI_Theta*ncol;
	for (uint coli = 0; coli < ncol; ++coli)
		{
		uint posQi = posQs[coli];
		uint posTi = posTs[coli];
		assert(posQi != UINT_MAX);
		assert(posTi != UINT_MAX);
		assert(posQi < LQ);
		assert(posTi < LT);

		for (uint colj = coli + 1; colj < ncol; ++colj)
			{
			uint posQj = posQs[colj];
			uint posTj = posTs[colj];
			assert(posQj != UINT_MAX);
			assert(posTj != UINT_MAX);
			assert(posQj < LQ);
			assert(posTj < LT);

			int diffij_Q = abs(int(posQi) - int(posQj));
			int diffij_T = abs(int(posTi) - int(posTj));
			if (diffij_Q > int(M) || diffij_T > int(M))
				continue;
			uint kQ = banded_ij_to_k(posQi, posQj);
			uint kT = banded_ij_to_k(posTi, posTj);
			sid_t sid_dQ_squared = distmxQ[kQ];
			sid_t sid_dT_squared = distmxT[kT];
			float dQ_squared = sid2dist2(sid_dQ_squared);
			float dT_squared = sid2dist2(sid_dT_squared);
			float dQ = sqrtf(dQ_squared);
			float dT = sqrtf(dT_squared);

			double DALI_dpscorefun(double a, double b);
			score += (float) DALI_dpscorefun(dQ, dT);
			}
		}
	return score;
	}

float flat_get_dali(
	const string &labelQ, const string &labelT,
	const string &path,
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const sid_t *distmxQ,
	const sid_t *distmxT)
	{
	const uint pathlen = uint(path.size());

// @@TODO std::vector here
	vector<uint32_t> posQs;
	vector<uint32_t> posTs;
	path2posvecs(labelQ, labelT, path, loQ, LQ, loT, LT, posQs, posTs);
	return flat_get_dali2(
		loQ, LQ, loT, LT,
		posQs, posTs,
		distmxQ, distmxT);
	}

float flat_get_dali3(
	const flat_aligner &fa,
	const sid_t *distmxQ,
	const sid_t *distmxT)
	{
	return flat_get_dali(
		fa.m_labelQ, fa.m_labelT,
		string(fa.m_path_buffer),
		fa.m_loQ, fa.m_LQ,
		fa.m_loT, fa.m_LT,
		distmxQ, distmxT);
	}


float flat_get_dalix2(
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const vector<uint> &posQs,
	const vector<uint> &posTs,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	float *colscores)
	{
	const uint ncol = uint(posQs.size());
	const uint M = flat_params::m_distmx_bandwidth;

	uint nlo = 0;
	for (uint coli = 0; coli < ncol; ++coli)
		{
		uint posQi = posQs[coli];
		uint posTi = posTs[coli];
		assert(posQi != UINT_MAX);
		assert(posTi != UINT_MAX);
		assert(posQi < LQ);
		assert(posTi < LT);

		float colscore = 0;
		for (uint colj = 0; colj < ncol; ++colj)
			{
			if (coli == colj)
				continue;
			uint posQj = posQs[colj];
			uint posTj = posTs[colj];
			assert(posQj != UINT_MAX);
			assert(posTj != UINT_MAX);
			assert(posQj < LQ);
			assert(posTj < LT);

			int diffij_Q = abs(int(posQi) - int(posQj));
			int diffij_T = abs(int(posTi) - int(posTj));
			if (diffij_Q > int(M) || diffij_T > int(M))
				continue;
			uint kQ = banded_ij_to_k(posQi, posQj);
			uint kT = banded_ij_to_k(posTi, posTj);
			sid_t sid_dQ_squared = distmxQ[kQ];
			sid_t sid_dT_squared = distmxT[kT];
			float dQ_squared = sid2dist2(sid_dQ_squared);
			float dT_squared = sid2dist2(sid_dT_squared);
			float dQ = sqrtf(dQ_squared);
			float dT = sqrtf(dT_squared);

			double DALI_dpscorefun(double a, double b);
			colscore += (float) DALI_dpscorefun(dQ, dT);
			}
		if (colscore <= tx)
			nlo += 1;
		colscores[coli] = colscore;
		}
	uint K = ncol - nlo;

	float score = 0; // s_DALI_Theta*K;
	for (uint coli = 0; coli < ncol; ++coli)
		{
		if (colscores[coli] <= tx)
			continue;
		uint posQi = posQs[coli];
		uint posTi = posTs[coli];
		assert(posQi != UINT_MAX);
		assert(posTi != UINT_MAX);
		assert(posQi < LQ);
		assert(posTi < LT);
		for (uint colj = coli + 1; colj < ncol; ++colj)
			{
			if (colscores[colj] <= tx)
				continue;
			uint posQj = posQs[colj];
			uint posTj = posTs[colj];
			assert(posQj != UINT_MAX);
			assert(posTj != UINT_MAX);
			assert(posQj < LQ);
			assert(posTj < LT);

			int diffij_Q = abs(int(posQi) - int(posQj));
			int diffij_T = abs(int(posTi) - int(posTj));
			if (diffij_Q > int(M) || diffij_T > int(M))
				continue;
			uint kQ = banded_ij_to_k(posQi, posQj);
			uint kT = banded_ij_to_k(posTi, posTj);
			sid_t sid_dQ_squared = distmxQ[kQ];
			sid_t sid_dT_squared = distmxT[kT];
			float dQ_squared = sid2dist2(sid_dQ_squared);
			float dT_squared = sid2dist2(sid_dT_squared);
			float dQ = sqrtf(dQ_squared);
			float dT = sqrtf(dT_squared);

			double DALI_dpscorefun(double a, double b);
			score += (float) DALI_dpscorefun(dQ, dT);
			}
		}

	return score;
	}

float flat_get_dalix(
	const char *path,
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	float *colscores)
	{
	Die("TODO");
//	const uint pathlen = (uint) strlen(path);
//
//// @@TODO std::vector here
//	vector<uint32_t> posQs;
//	vector<uint32_t> posTs;
//
//	path2posvecs(path, loQ, LQ, loT, LT, posQs, posTs);
//	return flat_get_dalix2(
//		loQ, LQ, loT, LT,
//		posQs, posTs,
//		distmxQ, distmxT, colscores);
	return 0;
	}

float flat_get_dalix3(
	const flat_aligner &fa,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	float *colscores)
	{
	Die("TODO");
	return 0;
	//return flat_get_dalix(
	//	fa.m_labelQ, fa.m_labelT,
	//	string(fa.m_path_buffer),
	//	fa.m_loQ, fa.m_LQ,
	//	fa.m_loT, fa.m_LT,
	//	distmxQ, distmxT,
	//	colscores);
	}

float flat_get_dalix(
	const flat_aligner &fa,
	const sid_t *distmxQ,
	const sid_t *distmxT)
	{
	Die("TODO");
	return 0;
	return flat_get_dali(
		fa.m_labelQ, fa.m_labelT,
		string(fa.m_path_buffer),
		fa.m_loQ, fa.m_loT,
		fa.m_LQ, fa.m_LT,
		distmxQ, distmxT);
	}
