#include "myutils.h"
#include "getticks.h"

#if DEBUG
static const uint SAMPLES = 10;
static const uint TIMING_SAMPLES = 10;
#else
static const uint SAMPLES = 1000;
static const uint TIMING_SAMPLES = 10000;
#endif

void read_profiles_and_logoddsmxvec(
	const string &specfn,
	vector<string> &feature_names,
	vector<uint> &alpha_sizes,
	vector<string> &labels,
	vector<vector<uint8_t> > &profiles,
	vector<vector<float> > &logoddsmxvec);

void check_profiles(
	vector<vector<uint8_t> > &profiles,
	vector<uint> &alpha_sizes);

uint32_t get_flat_pssm_feature_block_offsets(
	const uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	uint32_t * __restrict feature_block_offsets);

void fill_flat_pssm(
	const uint8_t * __restrict profQ,
	uint32_t LQ,
	uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	const uint32_t * __restrict feature_block_offsets,
	const float *const * __restrict weighted_logoddsmxvec,
	float * __restrict pssm);

void fill_smx_using_flat_pssm(
	const uint8_t * __restrict profA,
	uint32_t LA,
	uint32_t LB,
	uint32_t nfeat,
	const uint32_t * __restrict feature_block_offsets,
	const float * __restrict pssm,
	float * __restrict smx);

/***
Flattened PSSM.

For feature fi, a block of size alpha_sizes[fi] * LQ is stored:
    pssm_fi[a][j] = weighted_logoddsmxvec[fi][a*AS + profQ_fi[j]]

Layout is concatenation of feature blocks:
    [ feature 0 block ][ feature 1 block ] ... [ feature nfeat-1 block ]

Total number of floats required:
    LQ * sum_{fi=0}^{nfeat-1} alpha_sizes[fi]
***/

// Returns total number of floats required for flat PSSM.
// feature_block_offsets[fi] is the start offset, in floats, of feature fi.
uint32_t get_flat_pssm_feature_block_offsets(
	const uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	uint32_t * __restrict feature_block_offsets)
	{
	uint32_t offset = 0;
	for (uint32_t fi = 0; fi < nfeat; ++fi)
		{
		feature_block_offsets[fi] = offset;
		offset += alpha_sizes[fi];
		}
	return offset;
	}

/***
For feature fi, block starts at:
    pssm + feature_block_offsets[fi] * LQ

Within feature fi:
    pssm_fi[codeA][j] = weighted_logoddsmxvec[fi][codeA*AS + profQ_fi[j]]

feature_block_offsets are in units of rows, not floats.
Total floats required:
    LQ * sum_{fi=0}^{nfeat-1} alpha_sizes[fi]
***/
void fill_flat_pssm(
	const uint8_t * __restrict profQ,
	uint32_t LQ,
	uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	const uint32_t * __restrict feature_block_offsets,
	const float *const * __restrict weighted_logoddsmxvec,
	float * __restrict pssm)
	{
	for (uint32_t fi = 0; fi < nfeat; ++fi)
		{
		const uint32_t AS = alpha_sizes[fi];
		const float * __restrict mx = weighted_logoddsmxvec[fi];
		const uint8_t * __restrict profQ_fi = profQ + size_t(fi)*LQ;
		float * __restrict pssm_fi = pssm + size_t(feature_block_offsets[fi])*LQ;

		for (uint32_t a = 0; a < AS; ++a)
			{
			const float * __restrict mx_row = mx + size_t(a)*AS;
			float * __restrict pssm_row = pssm_fi + size_t(a)*LQ;

			for (uint32_t j = 0; j < LQ; ++j)
				pssm_row[j] = mx_row[profQ_fi[j]];
			}
		}
	}

/***
For feature fi, block starts at:
    pssm + feature_block_offsets[fi] * LB

Within feature fi:
    pssm_fi[codeA][j] is row codeA of length LB
***/
void fill_smx_using_flat_pssm(
	const uint8_t * __restrict profA,
	uint32_t LA,
	uint32_t LB,
	uint32_t nfeat,
	const uint32_t * __restrict feature_block_offsets,
	const float * __restrict pssm,
	float * __restrict smx)
	{
// Feature 0 handled as special case, note = instead of +=
	{
	const uint8_t * __restrict profA_f0 = profA;
	const float * __restrict pssm_f0 = pssm + size_t(feature_block_offsets[0])*LB;

	for (uint32_t i = 0; i < LA; ++i)
		{
		const uint8_t codeA_i = profA_f0[i];
		const float * __restrict pssm_row = pssm_f0 + size_t(codeA_i)*LB;
		float * __restrict smx_row = smx + size_t(i)*LB;

		for (uint32_t j = 0; j < LB; ++j)
			smx_row[j] = pssm_row[j];
		}
	}

	for (uint32_t fi = 1; fi < nfeat; ++fi)
		{
		const uint8_t * __restrict profA_fi = profA + size_t(fi)*LA;
		const float * __restrict pssm_fi = pssm + size_t(feature_block_offsets[fi])*LB;

		for (uint32_t i = 0; i < LA; ++i)
			{
			const uint8_t codeA_i = profA_fi[i];
			const float * __restrict pssm_row = pssm_fi + size_t(codeA_i)*LB;
			float * __restrict smx_row = smx + size_t(i)*LB;

			for (uint32_t j = 0; j < LB; ++j)
				smx_row[j] += pssm_row[j];
			}
		}
	}

void fill_smx(
	const uint8_t * __restrict profA,
	uint32_t LA,
	const uint8_t * __restrict profB,
	uint32_t LB,
	uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	const float *const * __restrict weighted_logoddsmxvec,
	float * __restrict smx)
	{

// Feature 0 handled as special case, note = instead of +=
	const uint32_t AS_0 = alpha_sizes[0];
	const float *weighted_logoddsmx_f0 = weighted_logoddsmxvec[0];
	const uint8_t *profA_f0 = profA;
	const uint8_t *profB_f0 = profB;
	for (uint32_t i = 0; i < LA; ++i)
		{
		const uint8_t codeA_i = profA_f0[i];
		const float *weighted_logoddsmx_row = weighted_logoddsmx_f0 + codeA_i*AS_0;
		float *smx_row = smx + i*LB;
		for (uint32_t j = 0; j < LB; ++j)
			{
			const uint8_t codeB_j = profB_f0[j];
			smx_row[j] = weighted_logoddsmx_row[codeB_j];
			}
		}

	for (uint32_t fi = 1; fi < nfeat; ++fi)
		{
		const uint32_t AS_fi = alpha_sizes[fi];
		const float *weighted_logoddsmx_fi = weighted_logoddsmxvec[fi];
		const uint8_t *profA_fi = profA + fi*LA;
		const uint8_t *profB_fi = profB + fi*LB;
		for (uint32_t i = 0; i < LA; ++i)
			{
			const uint8_t codeA_i = profA_fi[i];
			const float *weighted_logoddsmx_row = weighted_logoddsmx_fi + codeA_i*AS_fi;
			float *smx_row = smx + i*LB;
			for (uint32_t j = 0; j < LB; ++j)
				{
				const uint8_t codeB_j = profB_fi[j];
				smx_row[j] += weighted_logoddsmx_row[codeB_j];
				}
			}
		}
	}

static void fill_smx_slow(
	const vector<uint8_t> &profA,
	uint32_t LA,
	const vector<uint8_t> &profB,
	uint32_t LB,
	uint32_t nfeat,
	const vector<uint32_t> &alpha_sizes,
	const float *const *weighted_logoddsmxvec,
	float *smx)
	{
	asserta(SIZE(alpha_sizes) == nfeat);
	asserta(SIZE(profA) == nfeat*LA);
	asserta(SIZE(profB) == nfeat*LB);

	memset(smx, 0, LA*LB*sizeof(smx[0]));

	for (uint i = 0; i < LA; ++i)
		{
		for (uint j = 0; j < LB; ++j)
			{
			for (uint fi = 0; fi < nfeat; ++fi)
				{
				uint AS_fi = alpha_sizes[fi];
				uint8_t code_i = profA[fi*LA + i];
				uint8_t code_j = profB[fi*LB + j];
				asserta(code_i < AS_fi);
				asserta(code_j < AS_fi);
				smx[i*LB + j] += weighted_logoddsmxvec[fi][code_i*AS_fi + code_j];
				}
			}
		}
	}

static void cmp_smx(
	const float *smx1,
	const float *smx2,
	uint LA,
	uint LB)
	{
	for (uint i = 0; i < LA; ++i)
		for (uint j = 0; j < LB; ++j)
			if (!feq(smx1[i*LB + j], smx2[i*LB + j]))
				Die("smx[i=%u][j=%u] %.3g, %.3g",
					i, j, smx1[i*LB + j], smx2[i*LB + j]);
	}

void cmd_test_fill_smx()
	{
	const string &specfn = g_Arg1;
	vector<string> labels;
	vector<vector<uint8_t> > profiles;
	vector<string> feature_names;
	vector<uint> alpha_sizes;
	vector<vector<float> > logoddsmxvec;
	read_profiles_and_logoddsmxvec(
		specfn,
		feature_names,
		alpha_sizes,
		labels,
		profiles,
		logoddsmxvec);

	const uint nfeat = SIZE(feature_names);
	const uint nprof = SIZE(labels);

	asserta(SIZE(alpha_sizes) == nfeat);
	asserta(SIZE(profiles) == nprof);

	check_profiles(profiles, alpha_sizes);

	float **weighted_logoddsmxvec = myalloc(float *, nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = alpha_sizes[fi];
		asserta(AS >= 2 && AS < 256);
		weighted_logoddsmxvec[fi] = logoddsmxvec[fi].data();
		}

	uint32_t *feature_block_offsets = myalloc(uint32_t, nfeat);
	const uint32_t rows_per_pos =
		get_flat_pssm_feature_block_offsets(nfeat, alpha_sizes.data(), feature_block_offsets);

	for (uint sample = 0; sample < SAMPLES; ++sample)
		{
		uint idxA = randu32()%nprof;
		uint idxB = randu32()%nprof;

		const vector<uint8_t> &vector_profA = profiles[idxA];
		const vector<uint8_t> &vector_profB = profiles[idxB];

		uint LA = SIZE(vector_profA)/nfeat;
		uint LB = SIZE(vector_profB)/nfeat;

		const uint8_t *profA = vector_profA.data();
		const uint8_t *profB = vector_profB.data();

		float *smx = myalloc(float, LA*LB);
		float *smx_slow = myalloc(float, LA*LB);
		float *smx_pssm = myalloc(float, LA*LB);

		fill_smx(profA, LA, profB, LB, nfeat,
			alpha_sizes.data(), weighted_logoddsmxvec, smx);

		fill_smx_slow(vector_profA, LA, vector_profB, LB, nfeat,
			alpha_sizes, weighted_logoddsmxvec, smx_slow);

		cmp_smx(smx, smx_slow, LA, LB);
	
		uint nr_floats = LB*rows_per_pos;
		float *pssm = myalloc(float, nr_floats);

		fill_flat_pssm(profB, LB, nfeat, alpha_sizes.data(),
			feature_block_offsets, weighted_logoddsmxvec, pssm);

		fill_smx_using_flat_pssm(profA, LA, LB, nfeat,
			feature_block_offsets, pssm, smx_pssm);

		cmp_smx(smx, smx_pssm, LA, LB);
		}

	vector<uint> idxAs;
	vector<uint> idxBs;
	uint idxB = randu32()%nprof;
	for (uint sample = 0; sample < TIMING_SAMPLES; ++sample)
		{
		idxAs.push_back(randu32()%nprof);
		idxBs.push_back(idxB);
		if (sample%100 == 0)
			idxB = randu32()%nprof;
		}

	TICKS t1 = GetClockTicks();
	for (uint sample = 0; sample < TIMING_SAMPLES; ++sample)
		{
		uint idxA = idxAs[sample];
		uint idxB = idxBs[sample];

		const vector<uint8_t> &vector_profA = profiles[idxA];
		const vector<uint8_t> &vector_profB = profiles[idxB];

		uint LA = SIZE(vector_profA)/nfeat;
		uint LB = SIZE(vector_profB)/nfeat;

		const uint8_t *profA = vector_profA.data();
		const uint8_t *profB = vector_profB.data();

		float *smx_slow = myalloc(float, LA*LB);

		fill_smx_slow(vector_profA, LA, vector_profB, LB, nfeat,
			alpha_sizes, weighted_logoddsmxvec, smx_slow);
		}
	TICKS t2 = GetClockTicks();

	for (uint sample = 0; sample < TIMING_SAMPLES; ++sample)
		{
		uint idxA = idxAs[sample];
		uint idxB = idxBs[sample];

		const vector<uint8_t> &vector_profA = profiles[idxA];
		const vector<uint8_t> &vector_profB = profiles[idxB];

		uint LA = SIZE(vector_profA)/nfeat;
		uint LB = SIZE(vector_profB)/nfeat;

		const uint8_t *profA = vector_profA.data();
		const uint8_t *profB = vector_profB.data();

		float *smx = myalloc(float, LA*LB);

		fill_smx(profA, LA, profB, LB, nfeat,
			alpha_sizes.data(), weighted_logoddsmxvec, smx);
		}
	TICKS t3 = GetClockTicks();

	uint prev_idxB = UINT_MAX;
	float *pssm = 0;
	uint cached = 0;
	uint notcached = 0;
	for (uint sample = 0; sample < TIMING_SAMPLES; ++sample)
		{
		uint idxA = idxAs[sample];
		uint idxB = idxBs[sample];

		const vector<uint8_t> &vector_profA = profiles[idxA];
		const vector<uint8_t> &vector_profB = profiles[idxB];

		uint LA = SIZE(vector_profA)/nfeat;
		uint LB = SIZE(vector_profB)/nfeat;

		const uint8_t *profA = vector_profA.data();
		const uint8_t *profB = vector_profB.data();

		float *smx_pssm = myalloc(float, LA*LB);

		if (idxB == prev_idxB)
			++cached;
		else
			{
			++notcached;
			uint nr_floats = LB*rows_per_pos;
			pssm = myalloc(float, nr_floats);

			fill_flat_pssm(profB, LB, nfeat, alpha_sizes.data(),
				feature_block_offsets, weighted_logoddsmxvec, pssm);

			prev_idxB = idxB;
			}

		fill_smx_using_flat_pssm(profA, LA, LB, nfeat,
			feature_block_offsets, pssm, smx_pssm);
		}
	TICKS t4 = GetClockTicks();
	Progress("%u cached, %u not cached\n", cached, notcached);

	double tslow = double(t2 - t1);
	double tfast = double(t3 - t2);
	double tpssm = double(t4 - t3);

	ProgressLog("%s  slow\n", FloatToStr(tslow));
	ProgressLog("%s  fast\n", FloatToStr(tfast));
	ProgressLog("%s  pssm (%.1f%%)\n", FloatToStr(tpssm), GetPct(tfast-tpssm, tfast));
	}
