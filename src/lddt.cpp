#include "myutils.h"
#include "pdbchain.h"

static const float g_LDDT_R0 = 15;
static const float g_LDDT_R0_squared = g_LDDT_R0*g_LDDT_R0;
static const float g_LDDT_thresholds[4] = { 0.5, 1, 2, 4 };
static const uint g_nr_thresholds = 4;

double GetLDDT_mu(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  bool DaliScorerCompatible)
	{
	const uint nr_cols = SIZE(PosQs);
	if (nr_cols == 0)
		return 0;
	asserta(SIZE(PosTs) == nr_cols);
	float total = 0;
	uint nr_cols_considered = 0;
	for (uint coli = 0; coli < nr_cols; ++coli)
		{
		uint pos1i = PosQs[coli];
		uint pos2i = PosTs[coli];
		if (pos1i == UINT_MAX || pos2i == UINT_MAX)
			continue;

		++nr_cols_considered;
		uint nr_considered = 0;
		uint nr_preserved = 0;
		for (uint colj = 0; colj < nr_cols; ++colj)
			{
			if (coli == colj)
				continue;
			uint pos1j = PosQs[colj];
			uint pos2j = PosTs[colj];
			if (pos1j == UINT_MAX || pos2j == UINT_MAX)
				continue;

			float d1 = Q.GetDist(pos1i, pos1j);
			float d2 = T.GetDist(pos2i, pos2j);
			if (d1 > g_LDDT_R0 && d2 > g_LDDT_R0)
				continue;
			//if (coli == 912)
			//	Log("coli=%u colj=%u d1=%.1f d2=%.1f\n", coli, colj, d1, d2);//@@
			for (uint k = 0; k < g_nr_thresholds; ++k)
				{
				float t = g_LDDT_thresholds[k];
				nr_considered += 1;
				float diff = abs(d1 - d2);
				if (diff <= t)
					nr_preserved += 1;
				}
			}
		float score = 0;
		if (nr_considered > 0)
			score = float(nr_preserved)/nr_considered;
		total += score;
		//Log("coli %u preserved %u considered %u score %.4f\n",
		//	coli, nr_preserved, nr_considered, score);//@@
		}

	if (nr_cols_considered == 0)
		return 0;
	float avg = total/nr_cols_considered;
	return avg;
	}

double GetLDDT_mu_fast(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs)
	{
	const uint nr_cols = SIZE(PosQs);
	if (nr_cols == 0)
		return 0;
	asserta(SIZE(PosTs) == nr_cols);
	uint *nr_considered_vec = myalloc(uint, nr_cols);
	uint *nr_preserved_vec = myalloc(uint, nr_cols);
	zero_array(nr_considered_vec, nr_cols);
	zero_array(nr_preserved_vec, nr_cols);
	for (uint coli = 0; coli < nr_cols; ++coli)
		{
		uint pos1i = PosQs[coli];
		uint pos2i = PosTs[coli];
		assert(pos1i != UINT_MAX);
		assert(pos2i != UINT_MAX);

		for (uint colj = coli + 1; colj < nr_cols; ++colj)
			{
			uint pos1j = PosQs[colj];
			uint pos2j = PosTs[colj];
			assert(pos1j != UINT_MAX);
			assert(pos2j != UINT_MAX);

			float d1_squared = Q.GetDist2(pos1i, pos1j);
			float d2_squared = T.GetDist2(pos2i, pos2j);
			if (d1_squared > g_LDDT_R0_squared && d2_squared > g_LDDT_R0_squared)
				continue;

			float d1 = sqrtf(d1_squared);
			float d2 = sqrtf(d2_squared);
			for (uint k = 0; k < g_nr_thresholds; ++k)
				{
				float t = g_LDDT_thresholds[k];
				float diff = abs(d1 - d2);
				if (diff <= t)
					{
					nr_preserved_vec[coli] += 1;
					nr_preserved_vec[colj] += 1;
					}
				}
			nr_considered_vec[coli] += g_nr_thresholds;
			nr_considered_vec[colj] += g_nr_thresholds;
			}
		}

	float total = 0;
	for (uint col = 0; col < nr_cols; ++col)
		{
		float score = 0;
		uint nr_preserved = nr_preserved_vec[col];
		uint nr_considered = nr_considered_vec[col];
		if (nr_considered > 0)
			score = float(nr_preserved)/nr_considered;
		total += score;
		}
	myfree(nr_considered_vec);
	myfree(nr_preserved_vec);
	float avg = total/nr_cols;
	return avg;
	}
