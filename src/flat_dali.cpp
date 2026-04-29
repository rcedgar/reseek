#include "myutils.h"
#include "flat_chain.h"
#include "flat_distmx.h"
#include "flat_aligner.h"
#include "sort.h"

static const float tx = 1.25;//@@TODO param
//static float s_DALI_Theta = 0;

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
	const string &path,
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const sid_t *distmxQ,
	const sid_t *distmxT)
	{
	const uint pathlen = uint(path.size());

	vector<uint32_t> posQs;
	vector<uint32_t> posTs;

	void path2posvecs(
		const string &path,
		uint loQ, uint LQ,
		uint loT, uint LT,
		vector<uint> &posQs,
		vector<uint> &posTs);
	path2posvecs(path, loQ, LQ, loT, LT, posQs, posTs);
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
	return flat_get_dali(string(fa.m_path_buffer),
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
	const string &path,
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	float *colscores)
	{
	const uint pathlen = uint(path.size());

	vector<uint32_t> posQs;
	vector<uint32_t> posTs;

	void path2posvecs(
		const string &path,
		uint loQ, uint LQ,
		uint loT, uint LT,
		vector<uint> &posQs,
		vector<uint> &posTs);
	path2posvecs(path, loQ, LQ, loT, LT, posQs, posTs);
	return flat_get_dalix2(
		loQ, LQ, loT, LT,
		posQs, posTs,
		distmxQ, distmxT, colscores);
	}

float flat_get_dalix3(
	const flat_aligner &fa,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	float *colscores)
	{
	return flat_get_dalix(
		string(fa.m_path_buffer),
		fa.m_loQ, fa.m_LQ,
		fa.m_loT, fa.m_LT,
		distmxQ, distmxT,
		colscores);
	}

float flat_get_dalix(
	const flat_aligner &fa,
	const sid_t *distmxQ,
	const sid_t *distmxT)
	{
	return flat_get_dali(
		string(fa.m_path_buffer),
		fa.m_loQ, fa.m_loT,
		fa.m_LQ, fa.m_LT,
		distmxQ, distmxT);
	}
