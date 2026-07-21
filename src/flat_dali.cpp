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

static size_t dali_align_up(size_t n, size_t alignment)
	{
	asserta(alignment > 0);
	const size_t mask = alignment - 1;
	asserta((alignment & mask) == 0);
	return (n + mask) & ~mask;
	}

static float dali_pair_term(
	uint posQi, uint posTi,
	uint posQj, uint posTj,
	const sid_t *distmxQ, const sid_t *distmxT)
	{
	const uint M = flat_params::m_distmx_bandwidth;
	const uint diffij_Q = (posQi > posQj) ? (posQi - posQj) : (posQj - posQi);
	const uint diffij_T = (posTi > posTj) ? (posTi - posTj) : (posTj - posTi);
	if (diffij_Q > M || diffij_T > M)
		return 0.0f;

	const uint kQ = banded_ij_to_k(posQi, posQj);
	const uint kT = banded_ij_to_k(posTi, posTj);
	const float dQ = sqrtf(sid2dist2(distmxQ[kQ]));
	const float dT = sqrtf(sid2dist2(distmxT[kT]));
	return (float) DALI_dpscorefun(dQ, dT);
	}

size_t dali_greedy_terms_bytes(uint nmatch)
	{
	if (nmatch == 0)
		return 0;
	asserta(nmatch <= SIZE_T_MAX / nmatch);
	const size_t n2 = (size_t) nmatch * (size_t) nmatch;
	asserta(n2 <= (SIZE_T_MAX - (alignof(float) - 1)) / sizeof(float));
	return (alignof(float) - 1) + n2 * sizeof(float);
	}

size_t dali_greedy_work_bytes(uint nmatch)
	{
	if (nmatch == 0)
		return 0;
	asserta(nmatch <= (SIZE_T_MAX - (alignof(float) - 1)) / sizeof(float));
	size_t bytes = (alignof(float) - 1) + (size_t) nmatch * sizeof(float);
	asserta(bytes <= SIZE_T_MAX - (size_t) nmatch);
	return bytes + (size_t) nmatch;
	}

float flat_get_dali4(
	const uint *posQs, uint LQ, 
	const uint *posTs, uint LT, uint nmatch,
	const sid_t *distmxQ, const sid_t *distmxT)
	{
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

			score += dali_pair_term(
				posQi, posTi, posQj, posTj, distmxQ, distmxT);
			}
		}
	return score;
	}

float dali_greedy(
	const uint *posQs, uint LQ,
	const uint *posTs, uint LT, uint nmatch,
	const sid_t *distmxQ, const sid_t *distmxT,
	void *terms_scratch, size_t terms_bytes,
	void *work_scratch, size_t work_bytes,
	uint *retained_cols, uint retained_cols_capacity,
	uint &nretained)
	{
	nretained = 0;
	if (nmatch == 0)
		return 0.0f;

	asserta(posQs != 0);
	asserta(posTs != 0);
	asserta(distmxQ != 0);
	asserta(distmxT != 0);
	asserta(terms_scratch != 0);
	asserta(work_scratch != 0);
	asserta(retained_cols != 0);
	asserta(retained_cols_capacity >= nmatch);

	const size_t need_terms = dali_greedy_terms_bytes(nmatch);
	const size_t need_work = dali_greedy_work_bytes(nmatch);
	asserta(terms_bytes >= need_terms);
	asserta(work_bytes >= need_work);

	char *tp = (char *) terms_scratch;
	char *tend = tp + terms_bytes;
	tp = (char *) dali_align_up((size_t) tp, alignof(float));
	asserta((size_t) (tend - tp) >= (size_t) nmatch * (size_t) nmatch * sizeof(float));
	float *terms = (float *) tp;

	char *wp = (char *) work_scratch;
	char *wend = wp + work_bytes;
	wp = (char *) dali_align_up((size_t) wp, alignof(float));
	asserta((size_t) (wend - wp) >= (size_t) nmatch * sizeof(float));
	float *colscores = (float *) wp;
	wp += (size_t) nmatch * sizeof(float);
	asserta((size_t) (wend - wp) >= (size_t) nmatch);
	uint8_t *active = (uint8_t *) wp;

	extern float g_DALI_Theta;
	for (uint i = 0; i < nmatch; ++i)
		{
		const uint posQi = posQs[i];
		const uint posTi = posTs[i];
		asserta(posQi != UINT_MAX);
		asserta(posTi != UINT_MAX);
		asserta(posQi < LQ);
		asserta(posTi < LT);

		terms[(size_t) i * nmatch + i] = 0.0f;
		float cscore = g_DALI_Theta;
		for (uint j = 0; j < nmatch; ++j)
			{
			if (j == i)
				continue;
			const uint posQj = posQs[j];
			const uint posTj = posTs[j];
			asserta(posQj != UINT_MAX);
			asserta(posTj != UINT_MAX);
			asserta(posQj < LQ);
			asserta(posTj < LT);

			const float t = dali_pair_term(
				posQi, posTi, posQj, posTj, distmxQ, distmxT);
			terms[(size_t) i * nmatch + j] = t;
			cscore += t;
			}
		colscores[i] = cscore;
		active[i] = 1;
		}

	uint nactive = nmatch;
	for (;;)
		{
		uint worst = UINT_MAX;
		float worst_score = 0.0f;
		for (uint i = 0; i < nmatch; ++i)
			{
			if (!active[i])
				continue;
			if (colscores[i] >= 0.0f)
				continue;
			if (worst == UINT_MAX || colscores[i] < worst_score)
				{
				worst = i;
				worst_score = colscores[i];
				}
			}
		if (worst == UINT_MAX)
			break;

		active[worst] = 0;
		--nactive;
		const float *row = terms + (size_t) worst * nmatch;
		for (uint j = 0; j < nmatch; ++j)
			{
			if (!active[j])
				continue;
			colscores[j] -= row[j];
			}
		}

	nretained = 0;
	for (uint i = 0; i < nmatch; ++i)
		{
		if (!active[i])
			continue;
		retained_cols[nretained++] = i;
		}
	asserta(nretained == nactive);

	float score = g_DALI_Theta * (float) nretained;
	for (uint ii = 0; ii < nretained; ++ii)
		{
		const uint i = retained_cols[ii];
		for (uint jj = ii + 1; jj < nretained; ++jj)
			{
			const uint j = retained_cols[jj];
			score += terms[(size_t) i * nmatch + j];
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
