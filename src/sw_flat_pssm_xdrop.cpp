// Smith–Waterman flat PSSM with X-drop (banded + reference + tests).

#include "myutils.h"
#include <limits>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <atomic>
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
		uint8_t t;
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
			// SM first: seed / local start may also carry MI/MD gap bits.
			if (t & TRACEBITS_SM)
				{
				Leni = Besti - i + 1;
				Lenj = Bestj - j + 1;
				path_buffer[ncol] = 0;
				reverse_path_buffer(path_buffer, ncol);
				return;
				}
			if (t & TRACEBITS_DM) State = 'D';
			else if (t & TRACEBITS_IM) State = 'I';
			else State = 'M';
			--i; --j;
			break;
		case 'D':
			if (i == 0)
				{
				// Undo edge emit; alignment ends at previous cell.
				--ncol;
				path_buffer[ncol] = 0;
				Leni = Besti - i;
				Lenj = Bestj - j;
				reverse_path_buffer(path_buffer, ncol);
				return;
				}
			t = TB[(i - 1) * LB + j];
			State = (t & TRACEBITS_MD) ? 'M' : 'D';
			--i;
			break;
		case 'I':
			if (j == 0)
				{
				--ncol;
				path_buffer[ncol] = 0;
				Leni = Besti - i;
				Lenj = Bestj - j;
				reverse_path_buffer(path_buffer, ncol);
				return;
				}
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
	bool seed_anchored,
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
		// Bwd DP uses reverse indices: i=0,j=0 scores residues (posQ-1,posT-1).
		// Free (non-anchored) mode still bands from absolute high corner for
		// compatibility with existing tests; seed-anchored starts at (0,0).
		if (seed_anchored)
			{
			seedQ = 0;
			seedT = 0;
			}
		else
			{
			seedQ = posQ - 1;
			seedT = posT - 1;
			}
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

			uint8_t TraceBits = 0;
			SavedM0 = M0;

			// MATCH
			{
			float xM;
			const bool at_seed = (i_abs == seedQ && j == seedT);
			if (seed_anchored && at_seed)
				{
				// Force path to start with a match at the seed (Mu XDropHSP).
				xM = 0.0f;
				TraceBits = TRACEBITS_SM;
				}
			else
				{
				xM = M0;
				if (Drow[j] > xM) { xM = Drow[j]; TraceBits = TRACEBITS_DM; }
				if (I0 > xM) { xM = I0; TraceBits = TRACEBITS_IM; }
				if (!seed_anchored && 0.0f >= xM)
					{ xM = 0.0f; TraceBits = TRACEBITS_SM; }
				}

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
	bool seed_anchored,
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
		banded, XdropDir::Fwd, seed_anchored,
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
	return sw_flat_pssm_xdrop_fwd_impl(true, false, false,
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
	return sw_flat_pssm_xdrop_fwd_impl(true, true, false,
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
	return sw_flat_pssm_xdrop_fwd_impl(false, false, false,
		scratch_rows, TB, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext, active, loQ, loT, path_buffer, ncol);
	}

static float sw_flat_pssm_xdrop_bwd_impl(
	bool scoreonly,
	bool seed_anchored,
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
		true, XdropDir::Bwd, seed_anchored,
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
	return sw_flat_pssm_xdrop_bwd_impl(false, false,
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
	return sw_flat_pssm_xdrop_bwd_impl(true, false,
		scratch_rows, 0, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext, loQ, loT, 0, ncol);
	}

// True if path from (loQ,loT) stays inside [0,LQ) x [0,LT) on every match.
static bool xdrop_path_spans_ok(uint loQ, uint loT,
	const char *path, uint ncol, uint LQ, uint LT)
	{
	if (ncol == 0)
		return false;
	if (loQ >= LQ || loT >= LT)
		return false;
	uint q = loQ;
	uint t = loT;
	uint nmatch = 0;
	for (uint c = 0; c < ncol; ++c)
		{
		const char ch = path[c];
		if (ch == 'M')
			{
			if (q >= LQ || t >= LT)
				return false;
			++nmatch;
			++q;
			++t;
			}
		else if (ch == 'D')
			{
			++q;
			if (q > LQ)
				return false;
			}
		else if (ch == 'I')
			{
			++t;
			if (t > LT)
				return false;
			}
		else
			return false;
		}
	return nmatch > 0;
	}

// Strip leading/trailing I/D (adjust lo); require start/end M + valid spans.
static bool xdrop_canonicalize_local_path(uint &loQ, uint &loT,
	char *path, uint &ncol, uint LQ, uint LT, bool &did_trim)
	{
	did_trim = false;
	if (ncol == 0)
		return false;
	uint start = 0;
	while (start < ncol && (path[start] == 'I' || path[start] == 'D'))
		{
		if (path[start] == 'D')
			++loQ;
		else
			++loT;
		++start;
		did_trim = true;
		}
	uint end = ncol;
	while (end > start && (path[end - 1] == 'I' || path[end - 1] == 'D'))
		{
		--end;
		did_trim = true;
		}
	if (end <= start)
		return false;
	if (path[start] != 'M' || path[end - 1] != 'M')
		return false;
	const uint new_ncol = end - start;
	if (start > 0)
		memmove(path, path + start, new_ncol);
	path[new_ncol] = 0;
	ncol = new_ncol;
	return xdrop_path_spans_ok(loQ, loT, path, ncol, LQ, LT);
	}

// Checked lo = hi+1 - span; false on unsigned underflow.
static bool xdrop_bwd_lo_ok(uint hi_plus_1, uint span, uint &lo_out)
	{
	if (span > hi_plus_1)
		return false;
	lo_out = hi_plus_1 - span;
	return true;
	}

static std::atomic<uint> g_xdrop_hsp_calls{0};
static std::atomic<uint> g_xdrop_hsp_ok_merged{0};
static std::atomic<uint> g_xdrop_hsp_ok_fwd_only{0};
static std::atomic<uint> g_xdrop_hsp_ok_bwd_only{0};
static std::atomic<uint> g_xdrop_hsp_fail_empty{0};
static std::atomic<uint> g_xdrop_hsp_fail_span{0};
static std::atomic<uint> g_xdrop_hsp_fail_bwd_lo{0};
static std::atomic<uint> g_xdrop_hsp_fail_not_local{0};
static std::atomic<uint> g_xdrop_hsp_fwd_no_abut{0};
static std::atomic<uint> g_xdrop_hsp_merge_fallback_fwd{0};
static std::atomic<uint> g_xdrop_hsp_trimmed{0};

void reset_sw_flat_pssm_xdrop_hsp_stats()
	{
	g_xdrop_hsp_calls = 0;
	g_xdrop_hsp_ok_merged = 0;
	g_xdrop_hsp_ok_fwd_only = 0;
	g_xdrop_hsp_ok_bwd_only = 0;
	g_xdrop_hsp_fail_empty = 0;
	g_xdrop_hsp_fail_span = 0;
	g_xdrop_hsp_fail_bwd_lo = 0;
	g_xdrop_hsp_fail_not_local = 0;
	g_xdrop_hsp_fwd_no_abut = 0;
	g_xdrop_hsp_merge_fallback_fwd = 0;
	g_xdrop_hsp_trimmed = 0;
	}

void log_sw_flat_pssm_xdrop_hsp_stats()
	{
	ProgressLog("%10u  XDropHSP calls\n", g_xdrop_hsp_calls.load());
	ProgressLog("%10u  XDropHSP ok merged\n", g_xdrop_hsp_ok_merged.load());
	ProgressLog("%10u  XDropHSP ok fwd-only\n", g_xdrop_hsp_ok_fwd_only.load());
	ProgressLog("%10u  XDropHSP ok bwd-only\n", g_xdrop_hsp_ok_bwd_only.load());
	ProgressLog("%10u  XDropHSP fail empty\n", g_xdrop_hsp_fail_empty.load());
	ProgressLog("%10u  XDropHSP fail span\n", g_xdrop_hsp_fail_span.load());
	ProgressLog("%10u  XDropHSP fail bwd lo\n", g_xdrop_hsp_fail_bwd_lo.load());
	ProgressLog("%10u  XDropHSP fail not local (no terminal M)\n",
		g_xdrop_hsp_fail_not_local.load());
	ProgressLog("%10u  XDropHSP fwd no seed abut\n", g_xdrop_hsp_fwd_no_abut.load());
	ProgressLog("%10u  XDropHSP merge fell back to fwd\n",
		g_xdrop_hsp_merge_fallback_fwd.load());
	ProgressLog("%10u  XDropHSP path trimmed leading/trailing gaps\n",
		g_xdrop_hsp_trimmed.load());
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
	++g_xdrop_hsp_calls;

	const uint path_span = 2 * max(LQ, LT) + 4;
	char *fwd_path = path_buffer + path_span;
	char *bwd_path = path_buffer + 2 * path_span;
	uint fwd_ncol = 0, bwd_ncol = 0;
	uint fwd_loQ = 0, fwd_loT = 0, bwd_loQ = 0, bwd_loT = 0;

	const float fwd_score = sw_flat_pssm_xdrop_fwd_impl(true, false, true,
		scratch_rows, TB_fwd, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext, 0,
		fwd_loQ, fwd_loT, fwd_path, fwd_ncol);

	const float bwd_score = sw_flat_pssm_xdrop_bwd_impl(false, true,
		scratch_rows, TB_bwd, scratch_ppsms,
		profQ, LQ, pssmT, LT, feature_block_offsets, nfeat,
		posQ, posT, X, Open, Ext,
		bwd_loQ, bwd_loT, bwd_path, bwd_ncol);
	(void) bwd_loQ;
	(void) bwd_loT;

	auto fail = [&]() -> float
		{
		ncol = 0;
		path_buffer[0] = 0;
		loQ = loT = 0;
		return 0.0f;
		};

	auto finish_ok = [&](float score, uint which) -> float
		{
		bool did_trim = false;
		if (!xdrop_canonicalize_local_path(loQ, loT, path_buffer, ncol, LQ, LT, did_trim))
			{
			++g_xdrop_hsp_fail_not_local;
			return fail();
			}
		if (did_trim)
			++g_xdrop_hsp_trimmed;
		if (which == 0) ++g_xdrop_hsp_ok_merged;
		else if (which == 1) ++g_xdrop_hsp_ok_fwd_only;
		else ++g_xdrop_hsp_ok_bwd_only;
		return score;
		};

	auto take_fwd = [&]() -> float
		{
		if (fwd_ncol == 0 || fwd_score <= 0.0f)
			{
			++g_xdrop_hsp_fail_empty;
			return fail();
			}
		if (!xdrop_path_spans_ok(fwd_loQ, fwd_loT, fwd_path, fwd_ncol, LQ, LT))
			{
			++g_xdrop_hsp_fail_span;
			return fail();
			}
		memcpy(path_buffer, fwd_path, fwd_ncol + 1);
		ncol = fwd_ncol;
		loQ = fwd_loQ;
		loT = fwd_loT;
		return finish_ok(fwd_score, 1);
		};

	if (fwd_score + bwd_score <= 0.0f || (fwd_ncol == 0 && bwd_ncol == 0))
		{
		++g_xdrop_hsp_fail_empty;
		return fail();
		}

	if (bwd_ncol > 1)
		reverse_path_buffer(bwd_path, bwd_ncol);

	const bool fwd_abuts = (fwd_ncol == 0) ||
		(fwd_loQ == posQ && fwd_loT == posT);

	if (fwd_ncol > 0 && !fwd_abuts)
		{
		++g_xdrop_hsp_fwd_no_abut;
		return take_fwd();
		}

	if (bwd_ncol == 0)
		return take_fwd();

	if (fwd_ncol == 0)
		{
		uint BwdM, BwdD, BwdI;
		GetPathCounts(bwd_path, bwd_ncol, BwdM, BwdD, BwdI);
		uint lo_q = 0, lo_t = 0;
		if (!xdrop_bwd_lo_ok(posQ, BwdM + BwdD, lo_q) ||
			!xdrop_bwd_lo_ok(posT, BwdM + BwdI, lo_t))
			{
			++g_xdrop_hsp_fail_bwd_lo;
			return fail();
			}
		if (!xdrop_path_spans_ok(lo_q, lo_t, bwd_path, bwd_ncol, LQ, LT))
			{
			++g_xdrop_hsp_fail_span;
			return fail();
			}
		memcpy(path_buffer, bwd_path, bwd_ncol + 1);
		ncol = bwd_ncol;
		loQ = lo_q;
		loT = lo_t;
		return finish_ok(bwd_score, 2);
		}

	{
	uint BwdM, BwdD, BwdI;
	GetPathCounts(bwd_path, bwd_ncol, BwdM, BwdD, BwdI);
	uint merge_loQ = 0, merge_loT = 0;
	if (!xdrop_bwd_lo_ok(posQ, BwdM + BwdD, merge_loQ) ||
		!xdrop_bwd_lo_ok(posT, BwdM + BwdI, merge_loT))
		{
		++g_xdrop_hsp_fail_bwd_lo;
		++g_xdrop_hsp_merge_fallback_fwd;
		return take_fwd();
		}

	uint merged_ncol = 0;
	uint m_loQ = 0, m_loT = 0;
	MergeFwdBwdPaths(posQ, posT, fwd_path, fwd_ncol,
		posQ - 1, posT - 1, bwd_path, bwd_ncol,
		m_loQ, m_loT, path_buffer, merged_ncol);

	if (m_loQ != merge_loQ || m_loT != merge_loT ||
		!xdrop_path_spans_ok(m_loQ, m_loT, path_buffer, merged_ncol, LQ, LT))
		{
		++g_xdrop_hsp_fail_span;
		++g_xdrop_hsp_merge_fallback_fwd;
		return take_fwd();
		}

	ncol = merged_ncol;
	loQ = m_loQ;
	loT = m_loT;
	return finish_ok(fwd_score + bwd_score, 0);
	}
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

	float *scratch_rows = myalloc(float, 2 * LT + 3);
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
	float *scratch_rows = myalloc(float, 2 * LT + 3);
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

	// sw_flat_pssm init uses Drow[0..LT]; needs scratch_rows[2*LT+3]
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

	float *scratch_rows = myalloc(float, 2 * maxL + 3);
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
