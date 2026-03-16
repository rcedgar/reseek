#include "myutils.h"
#include "getticks.h"
#include "xdpmem.h"
#include "flat_sw_pssm.h"

#define USE_TPL	0

/***
$src/2026-03-15_benchmark_flat_fill_smx/2026-03-15_benchmark_flat_fill_smx

[721eae2]	2025-10_reseek_tune 
[1d210de]	reseek v2.9.i86linux64 
[1d210de]	reseek v2.9.win64

==> test_smx_CYGWIN.log <==
Cmp ticks 9.5e+08, malloc 25.8M (2.7%)
3.1G  slow
522.4M  fast
121.2M  pssm (76.8%)

==> test_smx_WSL.log <==
Cmp ticks 1.72e+09, malloc 57.9M (3.4%)
2.6G  slow
404.8M  fast
104.1M  pssm (74.3%)
***/

#if DEBUG
static const uint CMP_SAMPLES = 10;
static const uint TIMING_SAMPLES = 10;
static const uint SW_SAMPLES = 100;
#else
static const uint CMP_SAMPLES = 1000;
static const uint TIMING_SAMPLES = 10000;
static const uint SW_SAMPLES = 10000;
#endif
static const uint MAXL = 1024;

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

	float *smx_timing = myalloc(float, MAXL*MAXL);

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

	TICKS tcmp = GetClockTicks();
	TICKS tmalloc = 0;
	for (uint sample = 0; sample < CMP_SAMPLES; ++sample)
		{
		uint idxA = randu32()%nprof;
		uint idxB = randu32()%nprof;

		const vector<uint8_t> &vector_profA = profiles[idxA];
		const vector<uint8_t> &vector_profB = profiles[idxB];

		uint LA = SIZE(vector_profA)/nfeat;
		uint LB = SIZE(vector_profB)/nfeat;
		if (LA > MAXL || LB > MAXL)
			continue;

		const uint8_t *profA = vector_profA.data();
		const uint8_t *profB = vector_profB.data();

		TICKS t = GetClockTicks();
		float *smx = myalloc(float, LA*LB);
		float *smx_slow = myalloc(float, LA*LB);
		float *smx_pssm = myalloc(float, LA*LB);
		tmalloc += GetClockTicks() - t;

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
	tcmp = GetClockTicks() - tcmp;
	ProgressLog("Cmp ticks %.3g, malloc %s (%.1f%%)\n",
		double(tcmp), FloatToStr(double(tmalloc)), GetPct(double(tmalloc), double(tcmp)));

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
		if (LA > MAXL || LB > MAXL)
			continue;

		const uint8_t *profA = vector_profA.data();
		const uint8_t *profB = vector_profB.data();

		fill_smx_slow(vector_profA, LA, vector_profB, LB, nfeat,
			alpha_sizes, weighted_logoddsmxvec, smx_timing);
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
		if (LA > MAXL || LB > MAXL)
			continue;

		const uint8_t *profA = vector_profA.data();
		const uint8_t *profB = vector_profB.data();

		fill_smx(profA, LA, profB, LB, nfeat,
			alpha_sizes.data(), weighted_logoddsmxvec, smx_timing);
		}
	TICKS t3 = GetClockTicks();

	uint cached = 0;
	uint notcached = 0;
	{ // PSSM
	uint prev_idxB = UINT_MAX;

	uint nr_floats = MAXL*rows_per_pos;
	float *pssm = myalloc(float, nr_floats);
	for (uint sample = 0; sample < TIMING_SAMPLES; ++sample)
		{
		uint idxA = idxAs[sample];
		uint idxB = idxBs[sample];

		const vector<uint8_t> &vector_profA = profiles[idxA];
		const vector<uint8_t> &vector_profB = profiles[idxB];

		uint LA = SIZE(vector_profA)/nfeat;
		uint LB = SIZE(vector_profB)/nfeat;
		if (LA > MAXL || LB > MAXL)
			continue;

		const uint8_t *profA = vector_profA.data();
		const uint8_t *profB = vector_profB.data();

		if (idxB == prev_idxB)
			++cached;
		else
			{
			++notcached;

			fill_flat_pssm(profB, LB, nfeat, alpha_sizes.data(),
				feature_block_offsets, weighted_logoddsmxvec, pssm);

			prev_idxB = idxB;
			}

		fill_smx_using_flat_pssm(profA, LA, LB, nfeat,
			feature_block_offsets, pssm, smx_timing);
		}
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

/////////////////////////////////////////////////////////
// Smith-Waterman
/////////////////////////////////////////////////////////

static void cvt_smx(float *smx, uint LA, uint LB, float **SMxData)
	{
	for (uint i = 0; i < LA; ++i)
		SMxData[i] = smx + i*LB;
	}

float SWFast(XDPMem &Mem, const float * const *SMxData, uint LA, uint LB,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path);

static float Open = -10;
static float Ext = -1;
#include <float.h>	// FLT_MAX

/***
Scalar Smith-Waterman local alignment with affine gaps, using flat PSSM.

PSSM layout:
	For feature fi, block starts at:
		pssm + feature_block_offsets[fi] * LB

	Within feature fi:
		pssm_fi[codeA][j] is row codeA of length LB

Substitution score at cell (i,j):
	sum over fi of pssm_fi[ profA_fi[i] ][ j ]

Gap convention:
	gap_open   > 0
	gap_extend > 0

Recurrence:
	E(i,j) = max( H(i,j-1) - gap_open, E(i,j-1) - gap_extend )
	F(i,j) = max( H(i-1,j) - gap_open, F(i-1,j) - gap_extend )
	H(i,j) = max( 0, H(i-1,j-1) + sub(i,j), E(i,j), F(i,j) )

No traceback, returns best local score only.

Scratch requirements:
	H_prev: LB + 1 floats
	F_col:  LB + 1 floats

Caller must allocate these.
***/

float smith_waterman_affine_flat_pssm(
	const uint8_t * __restrict profA,
	uint32_t LA,
	uint32_t LB,
	uint32_t nfeat,
	const uint32_t * __restrict feature_block_offsets,
	const float * __restrict pssm,
	float gap_open,
	float gap_extend,
	float * __restrict H_prev,
	float * __restrict F_col)
	{
	asserta(gap_open >= 0.0f);
	asserta(gap_extend >= 0.0f);

	// Use a large negative finite number rather than -inf.
	// Safe because H is clamped to >= 0.
	const float NEG_INF = -FLT_MAX/4;

	// Row 0 initialization.
	H_prev[0] = 0.0f;
	F_col[0] = NEG_INF;
	for (uint32_t j = 1; j <= LB; ++j)
		{
		H_prev[j] = 0.0f;
		F_col[j] = NEG_INF;
		}

	float best_score = 0.0f;

	for (uint32_t i = 0; i < LA; ++i)
		{
		// Select PSSM row for each feature at this i.
		// If nfeat is small/fixed in practice, compiler may keep some of this in registers.
		// To avoid allocation, just re-walk features in the inner loop setup below.

		float H_diag = 0.0f;	// H(i-1,j-1)
		float H_left = 0.0f;	// H(i,j-1)
		float E = NEG_INF;		// E(i,j-1) carried across row

		H_prev[0] = 0.0f;		// H(i,0) after row update
		F_col[0] = NEG_INF;

		for (uint32_t j = 1; j <= LB; ++j)
			{
			const float H_up = H_prev[j];	// old H(i-1,j)

			// Compute substitution score at (i, j-1) from flat PSSM.
			float sub = 0.0f;
			for (uint32_t fi = 0; fi < nfeat; ++fi)
				{
				const uint8_t codeA_i = profA[size_t(fi)*LA + i];
				const float * __restrict pssm_fi =
					pssm + size_t(feature_block_offsets[fi])*LB;
				const float * __restrict pssm_row =
					pssm_fi + size_t(codeA_i)*LB;
				sub += pssm_row[j - 1];
				}

			// Horizontal gap state.
			{
			const float e_open = H_left - gap_open;
			const float e_ext = E - gap_extend;
			E = (e_open > e_ext ? e_open : e_ext);
			}

			// Vertical gap state.
			{
			const float f_open = H_up - gap_open;
			const float f_ext = F_col[j] - gap_extend;
			F_col[j] = (f_open > f_ext ? f_open : f_ext);
			}

			// Match/mismatch candidate.
			float H = H_diag + sub;

			if (E > H)
				H = E;
			if (F_col[j] > H)
				H = F_col[j];
			if (H < 0.0f)
				H = 0.0f;

			H_diag = H_up;
			H_prev[j] = H;
			H_left = H;

			if (H > best_score)
				best_score = H;
			}
		}

	return best_score;
	}
#include <float.h>

/***
Aggressively optimized scalar Smith-Waterman local alignment with affine gaps,
using flat PSSM directly. Score only, no traceback.

PSSM layout:
    For feature fi, block starts at:
        pssm + feature_block_offsets[fi] * LB

    Within feature fi:
        pssm_fi[codeA][j] is row codeA of length LB

Substitution score at cell (i,j):
    sum_fi pssm_fi[ profA_fi[i] ][ j ]

Gap convention:
    gap_open   > 0
    gap_extend > 0

DP states:
    E = horizontal gap state for current row
    F[j] = vertical gap state for column j
    H[j] = previous-row H on input, current-row H on output

Scratch required:
    H:        LB floats   (previous row, overwritten in place with current row)
    F:        LB floats   (vertical gap state)
    row_ptrs: nfeat pointers to currently selected PSSM rows

Caller allocates all scratch.
***/

static inline float max2f(float a, float b)
{
    return (a > b ? a : b);
}

static inline float max4f(float a, float b, float c, float d)
{
    float m = (a > b ? a : b);
    float n = (c > d ? c : d);
    return (m > n ? m : n);
}

float smith_waterman_affine_flat_pssm_fast(
    const uint8_t * __restrict profA,
    uint32_t LA,
    uint32_t LB,
    uint32_t nfeat,
    const uint32_t * __restrict feature_block_offsets,
    const float * __restrict pssm,
    float gap_open,
    float gap_extend,
    float * __restrict H,
    float * __restrict F,
    const float ** __restrict row_ptrs)
{
    asserta(gap_open >= 0.0f);
    asserta(gap_extend >= 0.0f);

    const float NEG_INF = -FLT_MAX/4;

    for (uint32_t j = 0; j < LB; ++j)
    {
        H[j] = 0.0f;
        F[j] = NEG_INF;
    }

    float best_score = 0.0f;

    for (uint32_t i = 0; i < LA; ++i)
    {
        // Select one PSSM row per feature for this i.
        for (uint32_t fi = 0; fi < nfeat; ++fi)
        {
            const uint8_t codeA_i = profA[size_t(fi)*LA + i];
            const float * __restrict pssm_fi =
                pssm + size_t(feature_block_offsets[fi]) * LB;
            row_ptrs[fi] = pssm_fi + size_t(codeA_i) * LB;
        }

        float H_left = 0.0f;      // H(i,j-1)
        float H_diag = 0.0f;      // H(i-1,j-1)
        float E = NEG_INF;        // E(i,j-1)

        switch (nfeat)
        {
        case 1:
            {
            const float * __restrict r0 = row_ptrs[0];
            for (uint32_t j = 0; j < LB; ++j)
            {
                const float H_up = H[j];
                const float sub = r0[j];

                const float e_open = H_left - gap_open;
                const float e_ext  = E      - gap_extend;
                E = (e_open > e_ext ? e_open : e_ext);

                const float f_open = H_up - gap_open;
                const float f_ext  = F[j] - gap_extend;
                F[j] = (f_open > f_ext ? f_open : f_ext);

                float h = H_diag + sub;
                if (E > h)    h = E;
                if (F[j] > h) h = F[j];
                if (h < 0.0f) h = 0.0f;

                H_diag = H_up;
                H[j] = h;
                H_left = h;
                if (h > best_score) best_score = h;
            }
            break;
            }

        case 2:
            {
            const float * __restrict r0 = row_ptrs[0];
            const float * __restrict r1 = row_ptrs[1];
            for (uint32_t j = 0; j < LB; ++j)
            {
                const float H_up = H[j];
                const float sub = r0[j] + r1[j];

                const float e_open = H_left - gap_open;
                const float e_ext  = E      - gap_extend;
                E = (e_open > e_ext ? e_open : e_ext);

                const float f_open = H_up - gap_open;
                const float f_ext  = F[j] - gap_extend;
                F[j] = (f_open > f_ext ? f_open : f_ext);

                float h = H_diag + sub;
                if (E > h)    h = E;
                if (F[j] > h) h = F[j];
                if (h < 0.0f) h = 0.0f;

                H_diag = H_up;
                H[j] = h;
                H_left = h;
                if (h > best_score) best_score = h;
            }
            break;
            }

        case 3:
            {
            const float * __restrict r0 = row_ptrs[0];
            const float * __restrict r1 = row_ptrs[1];
            const float * __restrict r2 = row_ptrs[2];
            for (uint32_t j = 0; j < LB; ++j)
            {
                const float H_up = H[j];
                const float sub = r0[j] + r1[j] + r2[j];

                const float e_open = H_left - gap_open;
                const float e_ext  = E      - gap_extend;
                E = (e_open > e_ext ? e_open : e_ext);

                const float f_open = H_up - gap_open;
                const float f_ext  = F[j] - gap_extend;
                F[j] = (f_open > f_ext ? f_open : f_ext);

                float h = H_diag + sub;
                if (E > h)    h = E;
                if (F[j] > h) h = F[j];
                if (h < 0.0f) h = 0.0f;

                H_diag = H_up;
                H[j] = h;
                H_left = h;
                if (h > best_score) best_score = h;
            }
            break;
            }

        case 4:
            {
            const float * __restrict r0 = row_ptrs[0];
            const float * __restrict r1 = row_ptrs[1];
            const float * __restrict r2 = row_ptrs[2];
            const float * __restrict r3 = row_ptrs[3];
            for (uint32_t j = 0; j < LB; ++j)
            {
                const float H_up = H[j];
                const float sub = r0[j] + r1[j] + r2[j] + r3[j];

                const float e_open = H_left - gap_open;
                const float e_ext  = E      - gap_extend;
                E = (e_open > e_ext ? e_open : e_ext);

                const float f_open = H_up - gap_open;
                const float f_ext  = F[j] - gap_extend;
                F[j] = (f_open > f_ext ? f_open : f_ext);

                float h = H_diag + sub;
                if (E > h)    h = E;
                if (F[j] > h) h = F[j];
                if (h < 0.0f) h = 0.0f;

                H_diag = H_up;
                H[j] = h;
                H_left = h;
                if (h > best_score) best_score = h;
            }
            break;
            }

        case 5:
            {
            const float * __restrict r0 = row_ptrs[0];
            const float * __restrict r1 = row_ptrs[1];
            const float * __restrict r2 = row_ptrs[2];
            const float * __restrict r3 = row_ptrs[3];
            const float * __restrict r4 = row_ptrs[4];
            for (uint32_t j = 0; j < LB; ++j)
            {
                const float H_up = H[j];
                const float sub = r0[j] + r1[j] + r2[j] + r3[j] + r4[j];

                const float e_open = H_left - gap_open;
                const float e_ext  = E      - gap_extend;
                E = (e_open > e_ext ? e_open : e_ext);

                const float f_open = H_up - gap_open;
                const float f_ext  = F[j] - gap_extend;
                F[j] = (f_open > f_ext ? f_open : f_ext);

                float h = H_diag + sub;
                if (E > h)    h = E;
                if (F[j] > h) h = F[j];
                if (h < 0.0f) h = 0.0f;

                H_diag = H_up;
                H[j] = h;
                H_left = h;
                if (h > best_score) best_score = h;
            }
            break;
            }

        case 6:
            {
            const float * __restrict r0 = row_ptrs[0];
            const float * __restrict r1 = row_ptrs[1];
            const float * __restrict r2 = row_ptrs[2];
            const float * __restrict r3 = row_ptrs[3];
            const float * __restrict r4 = row_ptrs[4];
            const float * __restrict r5 = row_ptrs[5];
            for (uint32_t j = 0; j < LB; ++j)
            {
                const float H_up = H[j];
                const float sub = r0[j] + r1[j] + r2[j] + r3[j] + r4[j] + r5[j];

                const float e_open = H_left - gap_open;
                const float e_ext  = E      - gap_extend;
                E = (e_open > e_ext ? e_open : e_ext);

                const float f_open = H_up - gap_open;
                const float f_ext  = F[j] - gap_extend;
                F[j] = (f_open > f_ext ? f_open : f_ext);

                float h = H_diag + sub;
                if (E > h)    h = E;
                if (F[j] > h) h = F[j];
                if (h < 0.0f) h = 0.0f;

                H_diag = H_up;
                H[j] = h;
                H_left = h;
                if (h > best_score) best_score = h;
            }
            break;
            }

        case 7:
            {
            const float * __restrict r0 = row_ptrs[0];
            const float * __restrict r1 = row_ptrs[1];
            const float * __restrict r2 = row_ptrs[2];
            const float * __restrict r3 = row_ptrs[3];
            const float * __restrict r4 = row_ptrs[4];
            const float * __restrict r5 = row_ptrs[5];
            const float * __restrict r6 = row_ptrs[6];
            for (uint32_t j = 0; j < LB; ++j)
            {
                const float H_up = H[j];
                const float sub = r0[j] + r1[j] + r2[j] + r3[j] + r4[j] + r5[j] + r6[j];

                const float e_open = H_left - gap_open;
                const float e_ext  = E      - gap_extend;
                E = (e_open > e_ext ? e_open : e_ext);

                const float f_open = H_up - gap_open;
                const float f_ext  = F[j] - gap_extend;
                F[j] = (f_open > f_ext ? f_open : f_ext);

                float h = H_diag + sub;
                if (E > h)    h = E;
                if (F[j] > h) h = F[j];
                if (h < 0.0f) h = 0.0f;

                H_diag = H_up;
                H[j] = h;
                H_left = h;
                if (h > best_score) best_score = h;
            }
            break;
            }

        case 8:
            {
            const float * __restrict r0 = row_ptrs[0];
            const float * __restrict r1 = row_ptrs[1];
            const float * __restrict r2 = row_ptrs[2];
            const float * __restrict r3 = row_ptrs[3];
            const float * __restrict r4 = row_ptrs[4];
            const float * __restrict r5 = row_ptrs[5];
            const float * __restrict r6 = row_ptrs[6];
            const float * __restrict r7 = row_ptrs[7];
            for (uint32_t j = 0; j < LB; ++j)
            {
                const float H_up = H[j];
                const float sub = r0[j] + r1[j] + r2[j] + r3[j] +
                                  r4[j] + r5[j] + r6[j] + r7[j];

                const float e_open = H_left - gap_open;
                const float e_ext  = E      - gap_extend;
                E = (e_open > e_ext ? e_open : e_ext);

                const float f_open = H_up - gap_open;
                const float f_ext  = F[j] - gap_extend;
                F[j] = (f_open > f_ext ? f_open : f_ext);

                float h = H_diag + sub;
                if (E > h)    h = E;
                if (F[j] > h) h = F[j];
                if (h < 0.0f) h = 0.0f;

                H_diag = H_up;
                H[j] = h;
                H_left = h;
                if (h > best_score) best_score = h;
            }
            break;
            }

        default:
            {
            for (uint32_t j = 0; j < LB; ++j)
            {
                const float H_up = H[j];

                float sub = 0.0f;
                uint32_t fi = 0;
                for (; fi + 3 < nfeat; fi += 4)
                    sub += row_ptrs[fi + 0][j] + row_ptrs[fi + 1][j]
                         + row_ptrs[fi + 2][j] + row_ptrs[fi + 3][j];
                for (; fi < nfeat; ++fi)
                    sub += row_ptrs[fi][j];

                const float e_open = H_left - gap_open;
                const float e_ext  = E      - gap_extend;
                E = (e_open > e_ext ? e_open : e_ext);

                const float f_open = H_up - gap_open;
                const float f_ext  = F[j] - gap_extend;
                F[j] = (f_open > f_ext ? f_open : f_ext);

                float h = H_diag + sub;
                if (E > h)    h = E;
                if (F[j] > h) h = F[j];
                if (h < 0.0f) h = 0.0f;

                H_diag = H_up;
                H[j] = h;
                H_left = h;
                if (h > best_score) best_score = h;
            }
            break;
            }
        }
    }

    return best_score;
}
void cmd_test_flat_sw()
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

	const float gap_open = -Open;
	const float gap_ext = -Ext;

	float *smx_timing = myalloc(float, MAXL*MAXL);
	float **SMxData = myalloc(float *, MAXL);
	float *H_prev = myalloc(float, MAXL+1);
	float *F_col = myalloc(float, MAXL+1);
	const float **row_ptrs = myalloc(const float *, nfeat);
	float *H = myalloc(float, MAXL);
	float *F = myalloc(float, MAXL);

	XDPMem Mem;

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

	double told = 0;
	double tpssm = 0;
	double tfast = 0;
#if USE_TPL
	double ttpl = 0;
#endif
	for (uint sample = 0; sample < SW_SAMPLES; ++sample)
		{
		uint idxA = randu32()%nprof;
		uint idxB = randu32()%nprof;

		if (sample%2 == 0)
			idxA = idxB;

		const vector<uint8_t> &vector_profA = profiles[idxA];
		const vector<uint8_t> &vector_profB = profiles[idxB];

		uint LA = SIZE(vector_profA)/nfeat;
		uint LB = SIZE(vector_profB)/nfeat;
		if (LA > MAXL || LB > MAXL)
			continue;

		const uint8_t *profA = vector_profA.data();
		const uint8_t *profB = vector_profB.data();

		TICKS t = GetClockTicks();
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

		TICKS t1 = GetClockTicks();
		float score = smith_waterman_affine_flat_pssm(
			profA, LA, LB, nfeat, feature_block_offsets, pssm,
			gap_open, gap_ext, H_prev, F_col);
		TICKS t2 = GetClockTicks();
		tpssm += double(t2 - t1);

		cmp_smx(smx, smx_pssm, LA, LB);

#if USE_TPL
		float score_tpl = 0;
		TICKS t5 = GetClockTicks();
		switch (nfeat)
			{
		case 2:
			score_tpl = smith_waterman_affine_flat_pssm_fixed_nfeat<2>(
				profA, LA, LB, feature_block_offsets, pssm,
				gap_open, gap_ext, H, F);
			break;

		case 6:
			score_tpl = smith_waterman_affine_flat_pssm_fixed_nfeat<6>(
				profA, LA, LB, feature_block_offsets, pssm,
				gap_open, gap_ext, H, F);
			break;

		default:
			Die("tpl %u", nfeat);
			}
		TICKS t6 = GetClockTicks();
		ttpl += double(t6 - t5);
#endif USE_TPL

		TICKS t7 = GetClockTicks();
		float score_fast = smith_waterman_affine_flat_pssm_fast(
			profA, LA, LB, nfeat, feature_block_offsets, pssm,
			gap_open, gap_ext, H, F, row_ptrs);
		TICKS t8 = GetClockTicks();
		tfast += double(t8 - t7);

		cvt_smx(smx, LA, LB, SMxData);

		uint Loi, Loj, Leni, Lenj;
		string Path;
		TICKS t3 = GetClockTicks();
		float score_old = SWFast(Mem, SMxData, LA, LB, Open, Ext,
			Loi, Loj, Leni, Lenj, Path);
		TICKS t4 = GetClockTicks();
		told += double(t4 - t3);

		//Log("%5u  score_old %8.3g   score %8.3g   %-20.20s...(%u)\n",
		//	sample, score_old, score, Path.c_str(), SIZE(Path));
		if (!feq(score, score_fast))
			Die("score %.3g fast %.3g", score, score_fast);
		if (!feq(score, score_old))
			Die("score %.3g old %.3g", score, score_old);
		if (!feq(score_fast, score_old))
			Die("score_fast %.3g old %.3g", score_fast, score_old);
#if USE_TPL
		if (!feq(score_tpl, score_old))
			Die("score_tpl %.3g old %.3g", score_tpl, score_old);
#endif
		}

	ProgressLog("%s  old\n", FloatToStr(told));
	ProgressLog("%s  pssm (%.1f%%)\n", FloatToStr(tpssm), 100 - GetPct(told-tpssm, told));
	ProgressLog("%s  fast (%.1f%%)\n", FloatToStr(tfast), 100 - GetPct(told-tfast, told));
#if USE_TPL
	ProgressLog("%s  tpl  (%.1f%%)\n", FloatToStr(ttpl), 100 - GetPct(told-ttpl, told));
#endif
	}
