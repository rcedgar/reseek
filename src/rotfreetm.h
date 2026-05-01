#pragma once

#include <stdint.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <string.h>

typedef uint16_t ic_t;
typedef uint16_t sid_t;
typedef uint32_t uint;

// -------------------------------------------------------------------
// Compact alignment representation.
// A matched pair t is (A[t], B[t]) with both arrays length K.
// Must be strictly increasing in both sequences.
// -------------------------------------------------------------------
struct align_path_t
{
	const uint32_t *A = 0;
	const uint32_t *B = 0;
	uint K = 0;

	bool IsValid() const
	{
		if (A == 0 || B == 0)
			return false;
		for (uint t = 1; t < K; ++t)
			{
			if (!(A[t] > A[t-1]))
				return false;
			if (!(B[t] > B[t-1]))
				return false;
			}
		return true;
	}
};

// -------------------------------------------------------------------
// Result container.
// per_pair_score[t] is the TM-like contribution of aligned pair t.
// delta_Ang[t] is the anchor-profile RMS discrepancy in Angstroms.
// anchors_* contain indexes into the alignment path, i.e. values in [0..K).
// -------------------------------------------------------------------
struct rotfree_tm_result_t
{
	float score = 0.0f;
	float sum = 0.0f;
	uint Lnorm = 0;

	std::vector<float> per_pair_score;
	std::vector<float> delta_Ang;

	std::vector<uint16_t> anchor_path_indexes;
};

// -------------------------------------------------------------------
// Parameters.
// -------------------------------------------------------------------
struct rotfree_tm_params_t
{
	// Number of anchors chosen from the alignment.
	uint anchor_count = 16;

	// Minimum spacing between chosen anchors in alignment-path coordinates.
	// Helps spread anchors across the alignment.
	uint min_anchor_sep = 8;

	// Robust clipping threshold in Angstroms for one anchor distance diff.
	// rho(x) = min(x*x, clip_Ang^2)
	float clip_Ang = 8.0f;

	// TM-like scale parameter.
	float d0_Ang = 3.0f;

	// If true, omit aligned pairs which are themselves anchors from scoring.
	// Usually false is fine.
	bool exclude_anchors_from_score = false;

	// Normalization length.
	// If 0, use alignment length K.
	uint Lnorm = 0;
};

// -------------------------------------------------------------------
// Utilities for sid_t.
//
// sid = (d_ic*d_ic)/16
// d_ic is in tenths of Angstroms
// d_Ang = 0.1*d_ic = sqrt(16*sid/100)
//       = 0.4*sqrt(sid)
//
// So squared distance in Angstrom^2 is:
// d2_Ang = 16*sid/100 = 0.16*sid
//
// Important:
// For score computation we never need sqrt at the single-anchor level.
// If diff_Ang^2 = (d1_Ang - d2_Ang)^2, then we do need d_Ang, hence sqrt.
// -------------------------------------------------------------------
static inline float sid_to_ang(sid_t sid)
{
	return 0.4f*sqrtf((float) sid);
}

static inline float sid_to_ang2(sid_t sid)
{
	return 0.16f*(float) sid;
}

static inline ic_t sid_to_ic(sid_t sid)
{
	float d_ic = 10.0f*sid_to_ang(sid);
	uint x = (uint) (d_ic + 0.5f);
	if (x > 65535u)
		x = 65535u;
	return (ic_t) x;
}

// -------------------------------------------------------------------
// Flat distance-matrix accessor.
//
// Stored for pairs with 1 <= j-i <= M.
// Layout:
//   distmx[M*i + (j - i - 1)]
//
// i,j are 0-based residue indexes.
// Caller is responsible for ensuring j > i and j-i <= M.
// -------------------------------------------------------------------
struct flat_sidmx_t
{
	const sid_t *m_Data = 0;
	uint m_L = 0;
	uint m_M = 0;

	flat_sidmx_t() {}
	flat_sidmx_t(const sid_t *Data, uint L, uint M)
		: m_Data(Data), m_L(L), m_M(M)
	{
	}

	inline bool HasPair(uint i, uint j) const
	{
		if (i == j)
			return false;
		if (i > j)
			{
			uint tmp = i;
			i = j;
			j = tmp;
			}
		uint d = j - i;
		return d >= 1 && d <= m_M && j < m_L;
	}

	inline sid_t Get(uint i, uint j) const
	{
		asserta(i != j);
		if (i > j)
			{
			uint tmp = i;
			i = j;
			j = tmp;
			}
		asserta(j < m_L);
		asserta(j > i);
		asserta(j - i <= m_M);
		return m_Data[m_M*i + (j - i - 1)];
	}
};

// -------------------------------------------------------------------
// Anchor-based rotation-free TM-like scorer.
//
// Concept:
// For an aligned pair t = (A[t], B[t]), compare its distances to a set
// of anchor pairs r = (A[a_r], B[a_r]).
//
// delta_t^2 = average_r min( (dA - dB)^2, clip^2 )
// score_t   = 1 / (1 + delta_t^2 / d0^2)
//
// Final score = sum_t score_t / Lnorm
//
// This is fast, invariant to rigid motion, and naturally compatible
// with local refinement.
// -------------------------------------------------------------------
class rotfree_tm_scorer
{
public:
	rotfree_tm_scorer() {}

	// ---------------------------------------------------------------
	// Simplest interface.
	// Picks anchors automatically from the alignment path.
	// ---------------------------------------------------------------
	float ScoreAlignment(
		const flat_sidmx_t &DM_A,
		const flat_sidmx_t &DM_B,
		const align_path_t &Path,
		const rotfree_tm_params_t &Params) const
	{
		rotfree_tm_result_t Result;
		ScoreAlignment(DM_A, DM_B, Path, Params, Result);
		return Result.score;
	}

	// ---------------------------------------------------------------
	// Full interface.
	// ---------------------------------------------------------------
	void ScoreAlignment(
		const flat_sidmx_t &DM_A,
		const flat_sidmx_t &DM_B,
		const align_path_t &Path,
		const rotfree_tm_params_t &Params,
		rotfree_tm_result_t &Result) const
	{
		asserta(Path.IsValid());
		asserta(DM_A.m_Data != 0);
		asserta(DM_B.m_Data != 0);

		Result.score = 0.0f;
		Result.sum = 0.0f;
		Result.Lnorm = (Params.Lnorm == 0 ? Path.K : Params.Lnorm);
		Result.anchor_path_indexes.clear();

		if (Path.K == 0 || Result.Lnorm == 0)
			return;

		ChooseAnchors(Path, Params, Result.anchor_path_indexes);

		Result.per_pair_score.assign(Path.K, 0.0f);
		Result.delta_Ang.assign(Path.K, 0.0f);

		const float clip2 = Params.clip_Ang*Params.clip_Ang;
		const float d02 = Params.d0_Ang*Params.d0_Ang;

		float Sum = 0.0f;
		uint Counted = 0;

		for (uint t = 0; t < Path.K; ++t)
			{
			if (Params.exclude_anchors_from_score && IsAnchor(t, Result.anchor_path_indexes))
				continue;

			const uint i = Path.A[t];
			const uint j = Path.B[t];

			float sum_rho = 0.0f;
			uint n = 0;

			for (uint q = 0; q < Result.anchor_path_indexes.size(); ++q)
				{
				const uint a = Result.anchor_path_indexes[q];
				if (a == t)
					continue;

				const uint ia = Path.A[a];
				const uint ja = Path.B[a];

				if (!DM_A.HasPair(i, ia))
					continue;
				if (!DM_B.HasPair(j, ja))
					continue;

				const sid_t sidA = DM_A.Get(i, ia);
				const sid_t sidB = DM_B.Get(j, ja);

				const float dA = sid_to_ang(sidA);
				const float dB = sid_to_ang(sidB);
				const float diff = dA - dB;
				float e2 = diff*diff;
				if (e2 > clip2)
					e2 = clip2;

				sum_rho += e2;
				++n;
				}

			float delta2 = clip2;
			if (n > 0)
				delta2 = sum_rho / (float) n;

			const float s = 1.0f / (1.0f + delta2 / d02);

			Result.per_pair_score[t] = s;
			Result.delta_Ang[t] = sqrtf(delta2);

			Sum += s;
			++Counted;
			}

		Result.sum = Sum;
		Result.score = Sum / (float) Result.Lnorm;
	}

	// ---------------------------------------------------------------
	// Optional interface when anchor path indexes are precomputed.
	// This is useful inside iterative refinement loops.
	// ---------------------------------------------------------------
	float ScoreAlignmentWithAnchors(
		const flat_sidmx_t &DM_A,
		const flat_sidmx_t &DM_B,
		const align_path_t &Path,
		const rotfree_tm_params_t &Params,
		const uint16_t *AnchorPathIndexes,
		uint AnchorCount) const
	{
		asserta(Path.IsValid());
		if (Path.K == 0)
			return 0.0f;

		const float clip2 = Params.clip_Ang*Params.clip_Ang;
		const float d02 = Params.d0_Ang*Params.d0_Ang;
		const uint Lnorm = (Params.Lnorm == 0 ? Path.K : Params.Lnorm);
		if (Lnorm == 0)
			return 0.0f;

		float Sum = 0.0f;

		for (uint t = 0; t < Path.K; ++t)
			{
			if (Params.exclude_anchors_from_score &&
				IsAnchor(t, AnchorPathIndexes, AnchorCount))
				continue;

			const uint i = Path.A[t];
			const uint j = Path.B[t];

			float sum_rho = 0.0f;
			uint n = 0;

			for (uint q = 0; q < AnchorCount; ++q)
				{
				const uint a = AnchorPathIndexes[q];
				if (a == t)
					continue;

				const uint ia = Path.A[a];
				const uint ja = Path.B[a];

				if (!DM_A.HasPair(i, ia))
					continue;
				if (!DM_B.HasPair(j, ja))
					continue;

				const float dA = sid_to_ang(DM_A.Get(i, ia));
				const float dB = sid_to_ang(DM_B.Get(j, ja));
				const float diff = dA - dB;

				float e2 = diff*diff;
				if (e2 > clip2)
					e2 = clip2;

				sum_rho += e2;
				++n;
				}

			float delta2 = clip2;
			if (n > 0)
				delta2 = sum_rho / (float) n;

			const float s = 1.0f / (1.0f + delta2 / d02);
			Sum += s;
			}

		return Sum / (float) Lnorm;
	}

private:
	static bool IsAnchor(uint t, const std::vector<uint16_t> &Anchors)
	{
		for (uint i = 0; i < Anchors.size(); ++i)
			if (Anchors[i] == t)
				return true;
		return false;
	}

	static bool IsAnchor(uint t, const uint16_t *Anchors, uint N)
	{
		for (uint i = 0; i < N; ++i)
			if (Anchors[i] == t)
				return true;
		return false;
	}

	// ---------------------------------------------------------------
	// Simple deterministic anchor selection.
	//
	// Strategy:
	//   - spread roughly uniformly over the alignment path
	//   - enforce min_anchor_sep if possible
	//   - skip very near ends when alignment is long enough
	//
	// This is intentionally simple for version 1.
	// Later you could replace this with:
	//   - best per-pair initial scores
	//   - local core residues
	//   - iterative reweighting
	// ---------------------------------------------------------------
	static void ChooseAnchors(
		const align_path_t &Path,
		const rotfree_tm_params_t &Params,
		std::vector<uint16_t> &Anchors)
	{
		Anchors.clear();

		if (Path.K == 0 || Params.anchor_count == 0)
			return;

		if (Path.K <= Params.anchor_count)
			{
			for (uint t = 0; t < Path.K; ++t)
				Anchors.push_back((uint16_t) t);
			return;
			}

		uint edge = 0;
		if (Path.K >= 32)
			edge = 2;

		const uint first = edge;
		const uint last = Path.K - 1 - edge;
		if (first > last)
			{
			for (uint t = 0; t < Path.K; ++t)
				Anchors.push_back((uint16_t) t);
			return;
			}

		for (uint a = 0; a < Params.anchor_count; ++a)
			{
			uint t;
			if (Params.anchor_count == 1)
				t = (first + last)/2;
			else
				t = first + (uint) (((uint64_t) a*(last - first) + (Params.anchor_count - 1)/2) /
									(Params.anchor_count - 1));

			if (!Anchors.empty())
				{
				uint prev = Anchors.back();
				if (t <= prev)
					t = prev + 1;
				}

			if (!Anchors.empty() && Params.min_anchor_sep > 1)
				{
				uint prev = Anchors.back();
				if (t < prev + Params.min_anchor_sep)
					t = prev + Params.min_anchor_sep;
				}

			if (t > last)
				break;

			Anchors.push_back((uint16_t) t);
			}

		// Fallback: if min spacing was too strict and we got too few,
		// fill in by uniform sampling without spacing.
		if (Anchors.size() == 0)
			{
			for (uint a = 0; a < Params.anchor_count; ++a)
				{
				uint t = first + (uint) (((uint64_t) a*(last - first) + (Params.anchor_count - 1)/2) /
										 (Params.anchor_count - 1));
				if (t > last)
					t = last;
				Anchors.push_back((uint16_t) t);
				}
			}
	}
};
