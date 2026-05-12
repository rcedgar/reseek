#include "myutils.h"
#include "flat_chain.h"
#include "flat_distmx.h"
#include "flat_aligner.h"
#include "flat_helpers.h"

// TODO float -> sid_t

//static const float R0 = 15;
//static const float g_LDDT_R0_squared = R0*R0;
////static const float g_LDDT_thresholds[4] = { 0.5, 1, 2, 4 };
//static const float g_LDDT_thresholds[] = { 1.4f };
//
//static const uint g_nr_thresholds =
//	sizeof(g_LDDT_thresholds)/sizeof(g_LDDT_thresholds[0]);

float flat_getlddt_muscle_some_floats(
	const uint32_t *posQs,
	const uint32_t LQ,
	const uint32_t *posTs,
	const uint32_t LT,
	const uint ncol,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	uint32_t *nr_considered_vec,
	uint32_t *nr_preserved_vec)
	{
	if (ncol == 0)
		return 0;
	const float R0 = flat_params::m_LDDT_R0;
	const float R02 = R0*R0;
	const uint nr_thresholds = flat_params::m_LDDT_nr_thresholds;
	const uint M = flat_params::m_distmx_bandwidth;
	const float *LDDT_thresholds = flat_params::m_LDDT_thresholds;
	zero_array(nr_considered_vec, ncol);
	zero_array(nr_preserved_vec, ncol);
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
			if (dQ > R0 && dT > R0)
				continue;

			for (uint k = 0; k < nr_thresholds; ++k)
				{
				// sid_t t = g_LDDT_thresholds2[k];
				float t = LDDT_thresholds[k];
				//uint16_t diff = abs(uint16(dQ_squared) - uint16(dT_squared));
				float diff = fabs(dQ - dT);
				if (diff <= t)
					{
					nr_preserved_vec[coli] += 1;
					nr_preserved_vec[colj] += 1;
					}
				}
			nr_considered_vec[coli] += nr_thresholds;
			nr_considered_vec[colj] += nr_thresholds;
			}
		}

	float total = 0;
	for (uint col = 0; col < ncol; ++col)
		{
		float score = 0;
		uint nr_preserved = nr_preserved_vec[col];
		uint nr_considered = nr_considered_vec[col];
		if (nr_considered > 0)
			score = float(nr_preserved)/nr_considered;
		total += score;
		}
	float avg = total/ncol;
	float lddt = avg;
	return lddt;
	}

float flat_getlddt_muscle_some_floats2(
	const sid_t *distmxQ,
	const sid_t *distmxT,
	uint LQ, uint LT,
	const vector<uint32_t> &posQs,
	const vector<uint32_t> &posTs,
	const uint M)
	{
	const uint ncol2 = uint(posQs.size());
	uint32_t *nr_considered_vec = myalloc(uint32_t, ncol2);
	uint32_t *nr_preserved_vec = myalloc(uint32_t, ncol2);
	float lddt = 
		flat_getlddt_muscle_some_floats(
			posQs.data(), LQ, posTs.data(), LT, ncol2,
			distmxQ, distmxT, nr_considered_vec, nr_preserved_vec);
	myfree(nr_considered_vec);
	myfree(nr_preserved_vec);
	return lddt;
	}

float flat_getlddt_muscle_some_floats3(
	const string &labelQ, const string &labelT,
	const string &path,
	uint32_t loQ, uint32_t LQ,
	uint32_t loT, uint32_t LT,
	const sid_t *distmxQ,
	const sid_t *distmxT)
	{
	const uint ncol = uint(path.size());

	vector<uint32_t> posQs;
	vector<uint32_t> posTs;
	path2posvecs(labelQ, labelT, path, loQ, LQ, loT, LT, posQs, posTs);

	const uint ncol2 = uint(posQs.size());
	asserta(posTs.size() == ncol2);
	uint32_t *nr_considered_vec = myalloc(uint32_t, ncol2);
	uint32_t *nr_preserved_vec = myalloc(uint32_t, ncol2);
	float lddt = 
		flat_getlddt_muscle_some_floats(
			posQs.data(), LQ, posTs.data(), LT, ncol2,
			distmxQ, distmxT, nr_considered_vec, nr_preserved_vec);
	myfree(nr_considered_vec);
	myfree(nr_preserved_vec);
	return lddt;
	}

float flat_getlddt_muscle_some_floats4(
	const flat_aligner &fa,
	const sid_t *distmxQ,
	const sid_t *distmxT)
	{
	return flat_getlddt_muscle_some_floats3(
		fa.m_labelQ, fa.m_labelT,
		string(fa.m_path_buffer),
		fa.m_loQ, fa.m_LQ,
		fa.m_loT, fa.m_LT,
		distmxQ, distmxT);
	}