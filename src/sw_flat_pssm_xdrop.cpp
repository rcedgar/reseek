// Smith–Waterman flat PSSM with X-drop (banded + reference + tests).

#include "myutils.h"
#include <limits>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include "sw_flat_pssm_xdrop.h"
#include "flat_helpers.h"

static const float MINUS_INFINITY = -9e9f;

// Heuristic xdrop vs full SW: within 5% of full or absolute 0.05.
static const double g_near_optimal_rel = 0.05;
static const double g_near_optimal_abs = 0.05;

struct flat_xdrop_heuristic_bucket
	{
	uint64_t attempted = 0;
	uint64_t optimal = 0;
	uint64_t near_optimal = 0;
	uint64_t sub_optimal = 0;
	uint64_t xdrop_gt_full = 0;
	uint64_t full_lo_before_seed = 0;
	};

struct flat_xdrop_test_stats
	{
	uint64_t banded_ref_attempted = 0;
	uint64_t banded_ref_mismatch = 0;
	uint64_t path_score_mismatch = 0;

	uint64_t x_huge_attempted = 0;
	uint64_t x_huge_mismatch = 0;

	uint64_t xdrop_vs_full_attempted = 0;
	uint64_t optimal = 0;
	uint64_t near_optimal = 0;
	uint64_t sub_optimal = 0;
	uint64_t xdrop_gt_full = 0;

	flat_xdrop_heuristic_bucket seed00;
	flat_xdrop_heuristic_bucket random_seed;

	uint64_t active_paths_attempted = 0;
	uint64_t active_path_score_mismatch = 0;
	uint64_t active_enum_inflated = 0;
	};

static flat_xdrop_test_stats g_flat_xdrop_stats;

static void report_heuristic_bucket(
	const char *label,
	const flat_xdrop_heuristic_bucket &b)
	{
	ProgressLog("  --- %s ---\n", label);
	ProgressLog("%10llu  attempted\n", (unsigned long long) b.attempted);
	ProgressLog("%10llu  optimal (xdrop == full)\n",
		(unsigned long long) b.optimal);
	ProgressLog("%10llu  near_optimal\n",
		(unsigned long long) b.near_optimal);
	ProgressLog("%10llu  sub_optimal\n",
		(unsigned long long) b.sub_optimal);
	ProgressLog("%10llu  xdrop_gt_full (bug)\n",
		(unsigned long long) b.xdrop_gt_full);
	ProgressLog("%10llu  full_SW_start_before_seed\n",
		(unsigned long long) b.full_lo_before_seed);
	}

static void report_flat_xdrop_test_stats()
	{
	const flat_xdrop_test_stats &s = g_flat_xdrop_stats;
	ProgressLog("\n=== flat_xdrop test summary ===\n");
	ProgressLog("%10llu  banded_vs_ref attempted\n",
		(unsigned long long) s.banded_ref_attempted);
	ProgressLog("%10llu  banded_vs_ref mismatch (Die)\n",
		(unsigned long long) s.banded_ref_mismatch);
	ProgressLog("%10llu  path_score mismatch (Die)\n",
		(unsigned long long) s.path_score_mismatch);
	ProgressLog("%10llu  X_huge vs full attempted\n",
		(unsigned long long) s.x_huge_attempted);
	ProgressLog("%10llu  X_huge mismatch (Die)\n",
		(unsigned long long) s.x_huge_mismatch);
	ProgressLog("%10llu  xdrop vs full SW attempted (total)\n",
		(unsigned long long) s.xdrop_vs_full_attempted);
	ProgressLog("%10llu  optimal (total)\n",
		(unsigned long long) s.optimal);
	ProgressLog("%10llu  near_optimal (total)\n",
		(unsigned long long) s.near_optimal);
	ProgressLog("%10llu  sub_optimal (total)\n",
		(unsigned long long) s.sub_optimal);
	ProgressLog("%10llu  xdrop_gt_full (total, bug)\n",
		(unsigned long long) s.xdrop_gt_full);
	report_heuristic_bucket("seed (0,0), random X", s.seed00);
	report_heuristic_bucket("random seed, random X", s.random_seed);
	ProgressLog("%10llu  active_paths attempted\n",
		(unsigned long long) s.active_paths_attempted);
	ProgressLog("%10llu  active path_score mismatch (Die)\n",
		(unsigned long long) s.active_path_score_mismatch);
	ProgressLog("%10llu  active enum > ref (informational)\n",
		(unsigned long long) s.active_enum_inflated);
	}

static void classify_xdrop_vs_full(
	float xdrop, float full,
	flat_xdrop_heuristic_bucket &bucket)
	{
	++g_flat_xdrop_stats.xdrop_vs_full_attempted;
	++bucket.attempted;
	if (feq(xdrop, full))
		{
		++g_flat_xdrop_stats.optimal;
		++bucket.optimal;
		return;
		}
	if (xdrop > full + float(g_near_optimal_abs))
		{
		++g_flat_xdrop_stats.xdrop_gt_full;
		++bucket.xdrop_gt_full;
		Die("xdrop %.4g > full %.4g", xdrop, full);
		}
	const double shortfall = double(full) - double(xdrop);
	if (shortfall <= g_near_optimal_abs)
		{
		++g_flat_xdrop_stats.near_optimal;
		++bucket.near_optimal;
		return;
		}
	if (full > 0.0f && shortfall / double(full) <= g_near_optimal_rel)
		{
		++g_flat_xdrop_stats.near_optimal;
		++bucket.near_optimal;
		return;
		}
	++g_flat_xdrop_stats.sub_optimal;
	++bucket.sub_optimal;
	}

static void test_sw_flat_pssm_xdrop_vs_full_trial(
	float *scratch_rows,
	uint8_t *TB,
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X,
	float open, float ext,
	flat_xdrop_heuristic_bucket &bucket)
	{
	char path_x[256], path_f[256];
	uint ncolx, ncolf, loQx, loTx, loQf, loTf;

	const float sx = sw_flat_pssm_xdrop_fwd(
		scratch_rows, TB, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, open, ext, loQx, loTx, path_x, ncolx);

	const float sf = sw_flat_pssm(
		scratch_rows, TB, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		open, ext, loQf, loTf, path_f, ncolf);

	if (loQf < posQ || loTf < posT)
		++bucket.full_lo_before_seed;

	classify_xdrop_vs_full(sx, sf, bucket);
	}

static void test_sw_flat_pssm_xdrop_vs_full(
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	float open, float ext,
	uint ntrial_seed00,
	uint ntrial_random_seed)
	{
	float *scratch_rows = myalloc(float, 2 * LT + 3);
	const float **scratch_ppsms = myalloc(const float *, nfeat);
	uint8_t *TB = myalloc(uint8_t, LQ * LT);

	for (uint trial = 0; trial < ntrial_seed00; ++trial)
		{
		const float X = float(1 + randu32() % 20);
		test_sw_flat_pssm_xdrop_vs_full_trial(
			scratch_rows, TB, scratch_ppsms,
			profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
			0, 0, X, open, ext, g_flat_xdrop_stats.seed00);
		}

	for (uint trial = 0; trial < ntrial_random_seed; ++trial)
		{
		const uint posQ = randu32() % max(1u, LQ);
		const uint posT = randu32() % max(1u, LT);
		const float X = float(1 + randu32() % 20);
		test_sw_flat_pssm_xdrop_vs_full_trial(
			scratch_rows, TB, scratch_ppsms,
			profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
			posQ, posT, X, open, ext, g_flat_xdrop_stats.random_seed);
		}

	myfree(scratch_rows);
	myfree(scratch_ppsms);
	myfree(TB);
	}

// ---------------------------------------------------------------------------
// PSSM column score at absolute (i, j)
// ---------------------------------------------------------------------------
static inline float pssm_at(
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint i, uint j)
	{
	float s = 0.0f;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		const uint8_t codeA = profQ[size_t(fi) * LQ + i];
		const float *pssm_fi = pssmT + size_t(feature_block_offsets[fi]) * LT;
		s += pssm_fi[size_t(codeA) * LT + j];
		}
	return s;
	}

static inline void select_pssm_rows(
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint i)
	{
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		const uint8_t codeA = profQ[size_t(fi) * LQ + i];
		const float *pssm_fi = pssmT + size_t(feature_block_offsets[fi]) * LT;
		scratch_ppsms[fi] = pssm_fi + size_t(codeA) * LT;
		}
	}

static inline float pssm_col_from_rows(
	const float **scratch_ppsms, uint nfeat, uint j)
	{
	float s = 0.0f;
	for (uint fi = 0; fi < nfeat; ++fi)
		s += scratch_ppsms[fi][j];
	return s;
	}

// Reverse pass: query index qr = posQ-1-i, target tc = posT-1-j, i,j in [0,posQ)x[0,posT)
static inline void select_pssm_rows_rev(
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	uint i, uint j)
	{
	const uint qr = posQ - 1 - i;
	const uint tc = posT - 1 - j;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		const uint8_t codeA = profQ[size_t(fi) * LQ + qr];
		const float *pssm_fi = pssmT + size_t(feature_block_offsets[fi]) * LT;
		scratch_ppsms[fi] = pssm_fi + size_t(codeA) * LT + tc;
		}
	}

static inline float pssm_col_rev(
	const float **scratch_ppsms, uint nfeat)
	{
	float s = 0.0f;
	for (uint fi = 0; fi < nfeat; ++fi)
		s += scratch_ppsms[fi][0];
	return s;
	}

// ---------------------------------------------------------------------------
// Traceback (duplicate of sw.cpp if you do not link sw.cpp)
// ---------------------------------------------------------------------------
static void reverse_path_buffer(char *path_buffer, uint ncol)
	{
	for (uint k = 0; k < ncol / 2; ++k)
		{
		char t = path_buffer[k];
		path_buffer[k] = path_buffer[ncol - 1 - k];
		path_buffer[ncol - 1 - k] = t;
		}
	path_buffer[ncol] = 0;
	}

static void traceback_flat_local(
	const uint8_t *TB, uint LA, uint LB,
	uint Besti, uint Bestj,
	uint &Leni, uint &Lenj,
	char *path_buffer, uint &ncol)
	{
	Leni = 0;
	Lenj = 0;
	uint i = Besti;
	uint j = Bestj;
	ncol = 0;
	char State = 'M';
	for (;;)
		{
		path_buffer[ncol++] = State;
		byte t;
		switch (State)
			{
		case 'M':
			if (i == 0 || j == 0)
				{
				Leni = Besti - i;
				Lenj = Bestj - j;
				path_buffer[ncol] = 0;
				reverse_path_buffer(path_buffer, ncol);
				return;
				}
			t = TB[(i - 1) * LB + (j - 1)];
			if (t & TRACEBITS_DM) State = 'D';
			else if (t & TRACEBITS_IM) State = 'I';
			else if (t & TRACEBITS_SM)
				{
				Leni = Besti - i + 1;
				Lenj = Bestj - j + 1;
				path_buffer[ncol] = 0;
				reverse_path_buffer(path_buffer, ncol);
				return;
				}
			else State = 'M';
			--i; --j;
			break;
		case 'D':
			asserta(i > 0);
			t = TB[(i - 1) * LB + j];
			State = (t & TRACEBITS_MD) ? 'M' : 'D';
			--i;
			break;
		case 'I':
			asserta(j > 0);
			t = TB[i * LB + (j - 1)];
			State = (t & TRACEBITS_MI) ? 'M' : 'I';
			--j;
			break;
		default:
			Die("traceback_flat_local: bad state %c", State);
			}
		}
	}

static void GetPathCounts(const char *Path, uint ncol, uint &M, uint &D, uint &I)
	{
	M = D = I = 0;
	for (uint c = 0; c < ncol; ++c)
		{
		switch (Path[c])
			{
		case 'M': ++M; break;
		case 'D': ++D; break;
		case 'I': ++I; break;
		default: Die("GetPathCounts: '%c'", Path[c]);
			}
		}
	}

static void MergeFwdBwdPaths(
	uint FwdLoQ, uint FwdLoT, const char *FwdPath, uint FwdNcol,
	uint BwdHiQ, uint BwdHiT, const char *BwdPath, uint BwdNcol,
	uint &loQ, uint &loT, char *out_path, uint &ncol)
	{
	asserta(FwdNcol > 0 || BwdNcol > 0);
	asserta(FwdLoQ == BwdHiQ + 1);
	asserta(FwdLoT == BwdHiT + 1);

	uint loA = FwdLoQ;
	uint loB = FwdLoT;
	if (BwdNcol > 0)
		{
		uint BwdM, BwdD, BwdI;
		GetPathCounts(BwdPath, BwdNcol, BwdM, BwdD, BwdI);
		loA = BwdHiQ + 1 - (BwdM + BwdD);
		loB = BwdHiT + 1 - (BwdM + BwdI);
		}
	ncol = 0;
	for (uint k = 0; k < BwdNcol; ++k)
		out_path[ncol++] = BwdPath[k];
	for (uint k = 0; k < FwdNcol; ++k)
		out_path[ncol++] = FwdPath[k];
	out_path[ncol] = 0;
	loQ = loA;
	loT = loB;
	}

// ---------------------------------------------------------------------------
// Core X-drop forward fill (banded or reference-unbanded)
// dir_fwd: true  => explore qr >= seedQ, tc >= seedT (global coords)
//          false => reverse extension, seed at (posQ-1,posT-1), qr = posQ-1-i
// ---------------------------------------------------------------------------
struct XdropFillResult
	{
	float best_score = 0.0f;
	uint best_i = UINT_MAX;
	uint best_j = UINT_MAX;
	};

enum class XdropDir { Fwd, Bwd };

static float sw_flat_pssm_xdrop_fill(
	bool banded,
	XdropDir dir,
	float *scratch_rows,
	uint8_t *TB,
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext,
	uint8_t *active,
	XdropFillResult &res)
	{
	res.best_score = 0.0f;
	res.best_i = UINT_MAX;
	res.best_j = UINT_MAX;

	asserta(Open <= 0.0f);
	asserta(Ext <= 0.0f);
	(void)banded;
	const float AbsOpen = -Open;
	const float AbsExt = -Ext;

	uint seedQ, seedT, i_min, i_max, j_min, j_max;
	if (dir == XdropDir::Fwd)
		{
		seedQ = posQ;
		seedT = posT;
		if (seedQ >= LQ || seedT >= LT)
			return 0.0f;
		i_min = seedQ;
		i_max = LQ - 1;
		j_min = seedT;
		j_max = LT - 1;
		}
	else
		{
		if (posQ == 0 || posT == 0)
			return 0.0f;
		seedQ = posQ - 1;
		seedT = posT - 1;
		i_min = 0;
		i_max = posQ - 1;
		j_min = 0;
		j_max = posT - 1;
		}

	float *Mrow = scratch_rows + 1;
	float *Drow = scratch_rows + LT + 2;
	Mrow[-1] = MINUS_INFINITY;

	for (uint j = 0; j < LT; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}

	float BestScore = 0.0f;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;

	uint prev_jlo = 0;
	uint prev_jhi = 0;
	uint jlo = seedT;
	uint jhi = seedT;

	float M0 = 0.0f;

	for (uint i_abs = i_min; i_abs <= i_max; ++i_abs)
		{
		if (dir == XdropDir::Fwd)
			select_pssm_rows(scratch_ppsms, profQ, LQ, pssmT, LT,
				feature_block_offsets, nfeat, i_abs);
		// Bwd: rows selected per (i_rev, j) inside inner loop

		if (jlo == prev_jlo)
			{
			if (jlo > 0)
				Mrow[jlo - 1] = MINUS_INFINITY;
			Drow[jlo] = MINUS_INFINITY;
			}

		uint endj = min(prev_jhi + 1, j_max);
		const uint band_extend_hi = min(jhi + 1, j_max);
		for (uint j = endj + 1; j <= band_extend_hi; ++j)
			{
			Mrow[j - 1] = MINUS_INFINITY;
			Drow[j] = MINUS_INFINITY;
			}

		uint next_jlo = UINT_MAX;
		uint next_jhi = UINT_MAX;
		float I0 = MINUS_INFINITY;

		uint8_t *TBrow = (TB != 0) ? (TB + size_t(i_abs) * LT) : 0;

		asserta(jlo <= jhi);

		float SavedM0 = MINUS_INFINITY;

		for (uint j = jlo; j <= jhi; ++j)
			{
			if (active)
				active[i_abs * LT + j] = 1;

			if (dir == XdropDir::Bwd)
				select_pssm_rows_rev(scratch_ppsms, profQ, LQ, pssmT, LT,
					feature_block_offsets, nfeat, posQ, posT,
					i_abs, j);

			byte TraceBits = 0;
			SavedM0 = M0;

			// MATCH
			{
			float xM = M0;
			if (Drow[j] > xM) { xM = Drow[j]; TraceBits = TRACEBITS_DM; }
			if (I0 > xM) { xM = I0; TraceBits = TRACEBITS_IM; }
			if (0.0f >= xM) { xM = 0.0f; TraceBits = TRACEBITS_SM; }

			M0 = Mrow[j];
			const float sub = (dir == XdropDir::Fwd)
				? pssm_col_from_rows(scratch_ppsms, nfeat, j)
				: pssm_col_rev(scratch_ppsms, nfeat);
			float s = xM + sub;
			Mrow[j] = s;

			float h = s - BestScore + X;
			if (h > 0.0f)
				{
				next_jlo = min(next_jlo, j + 1);
				next_jhi = j + 1;
				}
			if (h > AbsOpen)
				next_jlo = min(next_jlo, j);
			if (h > AbsExt && j == jhi && jhi + 1 <= j_max)
				{
				++jhi;
				uint new_endj = min(jhi + 1, j_max);
				new_endj = max(new_endj, endj);
				for (uint j2 = endj + 1; j2 <= new_endj; ++j2)
					{
					if (j2 - 1 > j)
						Mrow[j2 - 1] = MINUS_INFINITY;
					Drow[j2] = MINUS_INFINITY;
					}
				endj = new_endj;
				}

			if (s >= BestScore)
				{
				BestScore = s;
				Besti = i_abs;
				Bestj = j;
				}
			}

			// DELETE
			if (j != jlo)
				{
				float md = SavedM0 + Open;
				Drow[j] += Ext;
				if (md >= Drow[j])
					{
					Drow[j] = md;
					TraceBits |= TRACEBITS_MD;
					}
				float h = Drow[j] - BestScore + X;
				if (h > 0.0f)
					{
					next_jlo = min(next_jlo, j - 1);
					next_jhi = max(next_jhi, j - 1);
					}
				}

			// INSERT
			{
			float mi = SavedM0 + Open;
			I0 += Ext;
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				}
			float h = I0 - BestScore + X;
			if (h > 0.0f)
				{
				next_jlo = min(next_jlo, j + 1);
				next_jhi = max(next_jhi, j + 1);
				}
			if (h > AbsExt && j == jhi && jhi + 1 <= j_max)
				{
				++jhi;
				uint new_endj = min(jhi + 1, j_max);
				new_endj = max(new_endj, endj);
				for (uint j2 = endj + 1; j2 <= new_endj; ++j2)
					{
					Mrow[j2 - 1] = MINUS_INFINITY;
					Drow[j2] = MINUS_INFINITY;
					}
				endj = new_endj;
				}
			}

			if (TBrow)
				TBrow[j] = TraceBits;
			}

		// end Drow slot (XDropFwd special case)
		if (jhi < j_max)
			{
			const uint jhi1 = jhi + 1;
			if (TBrow)
				TBrow[jhi1] = 0;
			float md = M0 + Open;
			Drow[jhi1] += Ext;
			if (md >= Drow[jhi1])
				{
				Drow[jhi1] = md;
				if (TBrow)
					TBrow[jhi1] = TRACEBITS_MD;
				}
			}

		if (next_jlo == UINT_MAX)
			break;

		prev_jlo = jlo;
		prev_jhi = jhi;
		jlo = next_jlo;
		jhi = next_jhi;
		if (jlo > j_max)
			jlo = j_max;
		if (jhi > j_max)
			jhi = j_max;
		asserta(jlo <= jhi);
		asserta(jlo >= prev_jlo);

		if (jlo == prev_jlo)
			{
			M0 = MINUS_INFINITY;
			if (jlo <= j_max)
				Drow[jlo] = MINUS_INFINITY;
			}
		else
			M0 = Mrow[jlo - 1];
		}

	if (BestScore <= 0.0f)
		return 0.0f;

	res.best_score = BestScore;
	res.best_i = Besti;
	res.best_j = Bestj;
	return BestScore;
	}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------
static float sw_flat_pssm_xdrop_fwd_impl(
	bool banded,
	bool scoreonly,
	float *scratch_rows,
	uint8_t *TB,
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext,
	uint8_t *active,
	uint &loQ, uint &loT,
	char *path_buffer, uint &ncol)
	{
	loQ = loT = 0;
	ncol = 0;
	XdropFillResult res;
	const float score = sw_flat_pssm_xdrop_fill(
		banded, XdropDir::Fwd,
		scratch_rows, scoreonly ? 0 : TB,
		scratch_ppsms, profQ, LQ, pssmT, LT,
		feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext, active, res);
	if (score <= 0.0f || scoreonly)
		return score;

	uint Leni, Lenj;
	traceback_flat_local(TB, LQ, LT,
		res.best_i + 1, res.best_j + 1,
		Leni, Lenj, path_buffer, ncol);
	loQ = res.best_i + 1 - Leni;
	loT = res.best_j + 1 - Lenj;
	return score;
	}

float sw_flat_pssm_xdrop_fwd(
	float *scratch_rows, uint8_t *TB,
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext,
	uint &loQ, uint &loT,
	char *path_buffer, uint &ncol)
	{
	return sw_flat_pssm_xdrop_fwd_impl(true, false,
		scratch_rows, TB, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext, 0, loQ, loT, path_buffer, ncol);
	}

float sw_flat_pssm_xdrop_fwd_scoreonly(
	float *scratch_rows,
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext)
	{
	uint loQ, loT;
	uint ncol;
	return sw_flat_pssm_xdrop_fwd_impl(true, true,
		scratch_rows, 0, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext, 0, loQ, loT, 0, ncol);
	}

float sw_flat_pssm_xdrop_fwd_ref(
	float *scratch_rows, uint8_t *TB,
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext,
	uint8_t *active,
	uint &loQ, uint &loT,
	char *path_buffer, uint &ncol)
	{
	return sw_flat_pssm_xdrop_fwd_impl(false, false,
		scratch_rows, TB, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext, active, loQ, loT, path_buffer, ncol);
	}

static float sw_flat_pssm_xdrop_bwd_impl(
	bool scoreonly,
	float *scratch_rows, uint8_t *TB,
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext,
	uint &loQ, uint &loT,
	char *path_buffer, uint &ncol)
	{
	loQ = loT = 0;
	ncol = 0;
	if (posQ == 0 || posT == 0)
		return 0.0f;

	XdropFillResult res;
	const float score = sw_flat_pssm_xdrop_fill(
		true, XdropDir::Bwd,
		scratch_rows, scoreonly ? 0 : TB,
		scratch_ppsms, profQ, LQ, pssmT, LT,
		feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext, 0, res);
	if (score <= 0.0f || scoreonly)
		return score;

	uint Leni, Lenj;
	traceback_flat_local(TB, LQ, LT,
		res.best_i + 1, res.best_j + 1,
		Leni, Lenj, path_buffer, ncol);
	loQ = res.best_i + 1 - Leni;
	loT = res.best_j + 1 - Lenj;
	return score;
	}

float sw_flat_pssm_xdrop_bwd(
	float *scratch_rows, uint8_t *TB,
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext,
	uint &loQ, uint &loT,
	char *path_buffer, uint &ncol)
	{
	return sw_flat_pssm_xdrop_bwd_impl(false,
		scratch_rows, TB, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext, loQ, loT, path_buffer, ncol);
	}

float sw_flat_pssm_xdrop_bwd_scoreonly(
	float *scratch_rows,
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext)
	{
	uint loQ, loT;
	uint ncol;
	return sw_flat_pssm_xdrop_bwd_impl(true,
		scratch_rows, 0, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext, loQ, loT, 0, ncol);
	}

float sw_flat_pssm_xdrop_hsp(
	float *scratch_rows,
	uint8_t *TB_fwd, uint8_t *TB_bwd,
	const float **scratch_ppsms,
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext,
	uint &loQ, uint &loT,
	char *path_buffer, uint &ncol)
	{
	char *fwd_path = path_buffer;
	char *bwd_path = path_buffer + 2 * max(LQ, LT) + 4;
	uint fwd_ncol = 0, bwd_ncol = 0;
	uint fwd_loQ, fwd_loT, bwd_loQ, bwd_loT;

	const float fwd_score = sw_flat_pssm_xdrop_fwd(
		scratch_rows, TB_fwd, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext,
		fwd_loQ, fwd_loT, fwd_path, fwd_ncol);

	const float bwd_score = sw_flat_pssm_xdrop_bwd(
		scratch_rows, TB_bwd, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext,
		bwd_loQ, bwd_loT, bwd_path, bwd_ncol);

	if (fwd_score + bwd_score <= 0.0f)
		{
		ncol = 0;
		path_buffer[0] = 0;
		loQ = loT = 0;
		return 0.0f;
		}

	const uint FwdLoQ = (fwd_ncol > 0) ? fwd_loQ : posQ;
	const uint FwdLoT = (fwd_ncol > 0) ? fwd_loT : posT;
	const uint BwdHiQ = posQ - 1;
	const uint BwdHiT = posT - 1;
	MergeFwdBwdPaths(FwdLoQ, FwdLoT, fwd_path, fwd_ncol,
		BwdHiQ, BwdHiT, bwd_path, bwd_ncol,
		loQ, loT, path_buffer, ncol);
	return fwd_score + bwd_score;
	}

// ---------------------------------------------------------------------------
// Score a path (for validation; mirrors test_sw_enum.cpp)
// ---------------------------------------------------------------------------
static float score_path_pssm(
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint startQ, uint startT,
	float Open, float Ext,
	const char *path, uint ncol)
	{
	uint q = startQ, t = startT;
	float score = 0.0f;
	for (uint c = 0; c < ncol; ++c)
		{
		switch (path[c])
			{
		case 'M':
			score += pssm_at(nullptr, profQ, LQ, pssmT, LT,
				feature_block_offsets, nfeat, q, t);
			++q; ++t;
			break;
		case 'D':
			score += (c > 0 && path[c - 1] == 'M') ? Open : Ext;
			++q;
			break;
		case 'I':
			score += (c > 0 && path[c - 1] == 'M') ? Open : Ext;
			++t;
			break;
		default:
			Die("score_path_pssm");
			}
		}
	return score;
	}

// ---------------------------------------------------------------------------
// Tier-2: enumerate paths, keep those using only active M-cells, max score
// ---------------------------------------------------------------------------
static void enum_sw_paths_small(
	uint LA, uint LB,
	vector<uint> &starts_A,
	vector<uint> &starts_B,
	vector<string> &paths)
	{
	starts_A.clear();
	starts_B.clear();
	paths.clear();
	if (LA == 0 || LB == 0) return;

	string ops;
	function<void(uint, uint)> ExtendFromMatch;
	function<void(uint, uint)> ExtendFromGap;

	ExtendFromMatch = [&](uint a_used, uint b_used)
		{
		auto canonical = [](const string &s)
			{
			return s.find("DI") == string::npos && s.find("ID") == string::npos;
			};
		if (canonical(ops))
			{
			for (uint a0 = 0; a0 + a_used <= LA; ++a0)
				for (uint b0 = 0; b0 + b_used <= LB; ++b0)
					{
					starts_A.push_back(a0);
					starts_B.push_back(b0);
					paths.push_back(ops);
					}
			}
		if (a_used < LA && b_used < LB)
			{ ops.push_back('M'); ExtendFromMatch(a_used + 1, b_used + 1); ops.pop_back(); }
		if (a_used < LA)
			{ ops.push_back('D'); ExtendFromGap(a_used + 1, b_used); ops.pop_back(); }
		if (b_used < LB)
			{ ops.push_back('I'); ExtendFromGap(a_used, b_used + 1); ops.pop_back(); }
		};

	ExtendFromGap = [&](uint a_used, uint b_used)
		{
		if (a_used < LA && b_used < LB)
			{ ops.push_back('M'); ExtendFromMatch(a_used + 1, b_used + 1); ops.pop_back(); }
		if (a_used < LA)
			{ ops.push_back('D'); ExtendFromGap(a_used + 1, b_used); ops.pop_back(); }
		if (b_used < LB)
			{ ops.push_back('I'); ExtendFromGap(a_used, b_used + 1); ops.pop_back(); }
		};

	ExtendFromMatch(0, 0);
	}

static bool path_uses_only_active(
	const char *path, uint ncol,
	uint startQ, uint startT,
	const uint8_t *active, uint LT)
	{
	uint q = startQ, t = startT;
	for (uint c = 0; c < ncol; ++c)
		{
		if (path[c] == 'M')
			{
			if (!active[q * LT + t])
				return false;
			++q; ++t;
			}
		else if (path[c] == 'D')
			++q;
		else if (path[c] == 'I')
			++t;
		}
	return true;
	}

void test_sw_flat_pssm_xdrop_active_paths(
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float open, float ext)
	{
	if (LQ > 8 || LT > 8)
		return;

	++g_flat_xdrop_stats.active_paths_attempted;

	float *scratch_rows = myalloc(float, 2 * LT + 2);
	const float **scratch_ppsms = myalloc(const float *, nfeat);
	uint8_t *TB = myalloc(uint8_t, LQ * LT);
	uint8_t *active = myalloc(uint8_t, LQ * LT);
	memset(active, 0, LQ * LT);

	char path_buf[256];
	uint ncol, loQ, loT;

	const float ref_score = sw_flat_pssm_xdrop_fwd_ref(
		scratch_rows, TB, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, open, ext, active, loQ, loT, path_buf, ncol);

	if (ref_score > 0.0f)
		{
		const float path_sc = score_path_pssm(
			profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
			loQ, loT, open, ext, path_buf, ncol);
		if (!feq(ref_score, path_sc))
			{
			++g_flat_xdrop_stats.active_path_score_mismatch;
			Die("active_paths ref=%.4g path=%.4g seed=%u,%u X=%.2g",
				ref_score, path_sc, posQ, posT, X);
			}
		}

	vector<uint> starts_A, starts_B;
	vector<string> paths;
	enum_sw_paths_small(LQ, LT, starts_A, starts_B, paths);

	float best_enum = 0.0f;
	for (size_t p = 0; p < paths.size(); ++p)
		{
		if (starts_A[p] != loQ || starts_B[p] != loT)
			continue;
		const string &path = paths[p];
		if (!path_uses_only_active(path.c_str(), uint(path.size()),
			starts_A[p], starts_B[p], active, LT))
			continue;
		const float sc = score_path_pssm(
			profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
			starts_A[p], starts_B[p], open, ext,
			path.c_str(), uint(path.size()));
		if (sc > best_enum)
			best_enum = sc;
		}

	if (best_enum > ref_score + float(g_near_optimal_abs) &&
		!feq(best_enum, ref_score))
		++g_flat_xdrop_stats.active_enum_inflated;

	myfree(scratch_rows);
	myfree(scratch_ppsms);
	myfree(TB);
	myfree(active);
	}

void test_sw_flat_pssm_xdrop_banded_vs_ref(
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	float open, float ext,
	uint ntrial)
	{
	float *scratch_rows = myalloc(float, 2 * LT + 2);
	const float **scratch_ppsms = myalloc(const float *, nfeat);
	uint8_t *TBa = myalloc(uint8_t, LQ * LT);
	uint8_t *TBb = myalloc(uint8_t, LQ * LT);
	char path_a[256], path_b[256];

	for (uint trial = 0; trial < ntrial; ++trial)
		{
		++g_flat_xdrop_stats.banded_ref_attempted;

		const uint posQ = randu32() % max(1u, LQ);
		const uint posT = randu32() % max(1u, LT);
		const float X = float(1 + randu32() % 20);

		uint loQa, loTa, ncola, loQb, loTb, ncolb;
		memset(TBa, 0, LQ * LT);
		memset(TBb, 0, LQ * LT);

		const float sa = sw_flat_pssm_xdrop_fwd(
			scratch_rows, TBa, scratch_ppsms,
			profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
			posQ, posT, X, open, ext, loQa, loTa, path_a, ncola);

		const float sb = sw_flat_pssm_xdrop_fwd_ref(
			scratch_rows, TBb, scratch_ppsms,
			profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
			posQ, posT, X, open, ext, 0, loQb, loTb, path_b, ncolb);

		if (!feq(sa, sb))
			{
			++g_flat_xdrop_stats.banded_ref_mismatch;
			Die("banded vs ref score %.4g %.4g pos %u,%u X=%.2g",
				sa, sb, posQ, posT, X);
			}

		if (sa > 0.0f)
			{
			const float pa = score_path_pssm(
				profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
				loQa, loTa, open, ext, path_a, ncola);
			const float pb = score_path_pssm(
				profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
				loQb, loTb, open, ext, path_b, ncolb);
			if (!feq(sa, pa) || !feq(sb, pb))
				{
				++g_flat_xdrop_stats.path_score_mismatch;
				Die("path score mismatch fwd");
				}
			}
		}

	myfree(scratch_rows);
	myfree(scratch_ppsms);
	myfree(TBa);
	myfree(TBb);
	}

void test_sw_flat_pssm_xdrop_X_huge_vs_full(
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	float open, float ext)
	{
	++g_flat_xdrop_stats.x_huge_attempted;

	const float X = 1e9f;
	const uint posQ = 0, posT = 0;

	// sw_flat_pssm init uses Drow[0..LT]; needs scratch_rows[2*LT+2]
	float *scratch_rows = myalloc(float, 2 * LT + 3);
	const float **scratch_ppsms = myalloc(const float *, nfeat);
	uint8_t *TB = myalloc(uint8_t, LQ * LT);
	char path_x[256], path_f[256];
	uint ncolx, ncolf, loQx, loTx, loQf, loTf;

	const float sx = sw_flat_pssm_xdrop_fwd(
		scratch_rows, TB, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, open, ext, loQx, loTx, path_x, ncolx);

	const float sf = sw_flat_pssm(
		scratch_rows, TB, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		open, ext, loQf, loTf, path_f, ncolf);

	if (!feq(sx, sf))
		{
		++g_flat_xdrop_stats.x_huge_mismatch;
		Die("X_huge xdrop=%.4g full=%.4g", sx, sf);
		}

	myfree(scratch_rows);
	myfree(TB);
	}

void cmd_test_flat_xdrop()
	{
	memset(&g_flat_xdrop_stats, 0, sizeof(g_flat_xdrop_stats));

	const uint nfeat = 3;
	const uint minL = 3, maxL = 10, nprof = 100;
	uint32_t alpha_sizes[3] = { 3, 4, 5 };
	uint32_t feature_block_offsets[3];
	const uint sum_as = get_flat_pssm_feature_block_offsets(
		nfeat, alpha_sizes, feature_block_offsets);

	float *scratch_rows = myalloc(float, 2 * maxL + 2);
	const float **scratch_ppsms = myalloc(const float *, nfeat);
	uint8_t *TB = myalloc(uint8_t, maxL * maxL);
	float *pssm = myalloc(float, maxL * sum_as);

	float *wlo[nfeat];
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		wlo[fi] = myalloc(float, alpha_sizes[fi] * alpha_sizes[fi]);
		for (uint k = 0; k < alpha_sizes[fi] * alpha_sizes[fi]; ++k)
			wlo[fi][k] = float(int(randu32() % 7) - 3);
		}

	uint8_t *profs[nprof];
	uint32_t prof_len[nprof];
	for (uint i = 0; i < nprof; ++i)
		{
		prof_len[i] = minL + randu32() % (maxL - minL);
		profs[i] = myalloc(uint8_t, nfeat * prof_len[i]);
		for (uint fi = 0; fi < nfeat; ++fi)
			for (uint p = 0; p < prof_len[i]; ++p)
				profs[i][fi * prof_len[i] + p] =
					uint8_t(randu32() % alpha_sizes[fi]);
		}

	for (uint j = 0; j < nprof; ++j)
		{
		ProgressStep(j, nprof, "Testing");
		fill_flat_pssm(profs[j], prof_len[j], nfeat, alpha_sizes,
			feature_block_offsets, wlo, pssm);
		for (uint i = 0; i < nprof; ++i)
			{
			test_sw_flat_pssm_xdrop_banded_vs_ref(
				profs[i], prof_len[i], pssm, prof_len[j],
				feature_block_offsets, nfeat,
				-1.0f, -0.1f, 50);
			test_sw_flat_pssm_xdrop_X_huge_vs_full(
				profs[i], prof_len[i], pssm, prof_len[j],
				feature_block_offsets, nfeat, -1.0f, -0.1f);
			test_sw_flat_pssm_xdrop_vs_full(
				profs[i], prof_len[i], pssm, prof_len[j],
				feature_block_offsets, nfeat,
				-1.0f, -0.1f, 20, 20);
			for (uint pq = 0; pq < prof_len[i]; ++pq)
				for (uint pt = 0; pt < prof_len[j]; ++pt)
					test_sw_flat_pssm_xdrop_active_paths(
						profs[i], prof_len[i], pssm, prof_len[j],
						feature_block_offsets, nfeat,
						pq, pt, 5.0f, -1.0f, -0.1f);
			}
		}

	myfree(scratch_rows);
	myfree(TB);
	myfree(pssm);

	report_flat_xdrop_test_stats();
	}
