#pragma once

#include <stdint.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <string.h>

/***
Main idea

Anchors are a small set of aligned core residues that 
define a rotation-free geometric reference frame by 
comparing each candidate correspondence through its 
distances to those residues.

For aligned pair (i,j), instead of only comparing residue i in X to residue j in Y, allow the target in Y to be one of:

residue j
midpoint of segment (j−1, j)
midpoint of segment (j, j+1)

or more generally a small set of discrete segment positions.

The rotation-free discrepancy is still based on anchor-distance profiles.

For anchor pair a=(i_a, j_a) compare

d_X(i,i_a)

versus 

distance from the candidate point on Y backbone to anchor residue j_a

For a point on segment (j_0, ,j_1)
	​
at fraction λ in [0,1], distance squared to anchor k 
uses only inter-residue distances.
***/

#ifndef asserta
#include <assert.h>
#define asserta(x) assert(x)
#endif

typedef uint16_t ic_t;
typedef uint16_t sid_t;
typedef uint32_t uint;

// ---------------------------------------------------------------
// Alignment path: matched residue pairs.
// A[t], B[t], t=0..K-1
// ---------------------------------------------------------------
struct align_path_t
{
	const uint16_t* A = 0;
	const uint16_t* B = 0;
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

// ---------------------------------------------------------------
// Flat banded sid matrix.
// Stores pairs 1 <= |i-j| <= M
// layout: Data[M*i + (j-i-1)] for j=i+1..i+M
// ---------------------------------------------------------------
struct flat_sidmx_t
{
	const sid_t* m_Data = 0;
	uint m_L = 0;
	uint m_M = 0;

	flat_sidmx_t() {}
	flat_sidmx_t(const sid_t* Data, uint L, uint M)
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

// ---------------------------------------------------------------
// sid helpers
// d_Ang = sqrt(16*sid/100) = 0.4*sqrt(sid)
// d2_Ang = 0.16*sid
// ---------------------------------------------------------------
static inline float sid_to_ang(sid_t sid)
{
	return 0.4f*sqrtf((float) sid);
}

static inline float sid_to_ang2(sid_t sid)
{
	return 0.16f*(float) sid;
}

// ---------------------------------------------------------------
// Candidate target type in Y.
// ---------------------------------------------------------------
enum seg_target_type_t
{
	SEG_TARGET_RESIDUE = 0,   // residue j
	SEG_TARGET_SEGMENT = 1    // point on segment (j0,j1) at lambda_num/lambda_den
};

struct seg_target_t
{
	uint16_t kind = SEG_TARGET_RESIDUE;

	// for residue target: j0=j1=j
	// for segment target: segment endpoints j0,j1 with j1=j0+1 typically
	uint16_t j0 = 0;
	uint16_t j1 = 0;

	// lambda = lambda_num / lambda_den
	// residue target uses lambda = 0
	uint8_t lambda_num = 0;
	uint8_t lambda_den = 1;
};

// ---------------------------------------------------------------
// Result.
// chosen_target[t] tells which Y backbone point was best for pair t.
// delta_Ang[t] is RMS anchor-profile discrepancy.
// per_pair_score[t] is TM-like contribution.
// anchor_path_indexes are indexes into alignment path [0..K).
// ---------------------------------------------------------------
struct rotfree_tm_segment_result_t
{
	float score = 0.0f;
	float sum = 0.0f;
	uint Lnorm = 0;

	std::vector<float> per_pair_score;
	std::vector<float> delta_Ang;
	std::vector<seg_target_t> chosen_target;
	std::vector<uint16_t> anchor_path_indexes;
};

// ---------------------------------------------------------------
// Parameters.
// ---------------------------------------------------------------
struct rotfree_tm_segment_params_t
{
	uint anchor_count = 16;
	uint min_anchor_sep = 8;

	float clip_Ang = 8.0f;
	float d0_Ang = 3.0f;

	// if 0 use Path.K
	uint Lnorm = 0;

	// Allowed discrete target positions in Y for aligned pair (i,j).
	// Common simple settings:
	//   allow_residue = true
	//   allow_left_midpoint = true
	//   allow_right_midpoint = true
	bool allow_residue = true;
	bool allow_left_midpoint = true;   // midpoint of (j-1,j)
	bool allow_right_midpoint = true;  // midpoint of (j,j+1)

	// Optional extra fractions on adjacent segments.
	// Example: {1,3} with denom=4 means test 1/4 and 3/4.
	// Fractions are used on both left and right adjacent segments when enabled.
	std::vector<uint8_t> extra_lambda_nums;
	uint8_t extra_lambda_den = 1;

	bool exclude_anchors_from_score = false;
};

// ---------------------------------------------------------------
// Rotation-free point-to-backbone scorer.
// Scores X residues against candidate points on Y backbone.
// ---------------------------------------------------------------
class rotfree_tm_segment_scorer
{
public:
	float ScoreAlignment(
		const flat_sidmx_t& DM_X,
		const flat_sidmx_t& DM_Y,
		const align_path_t& Path,
		const rotfree_tm_segment_params_t& Params) const
	{
		rotfree_tm_segment_result_t Result;
		ScoreAlignment(DM_X, DM_Y, Path, Params, Result);
		return Result.score;
	}

	void ScoreAlignment(
		const flat_sidmx_t& DM_X,
		const flat_sidmx_t& DM_Y,
		const align_path_t& Path,
		const rotfree_tm_segment_params_t& Params,
		rotfree_tm_segment_result_t& Result) const
	{
		asserta(Path.IsValid());
		asserta(DM_X.m_Data != 0);
		asserta(DM_Y.m_Data != 0);

		Result.score = 0.0f;
		Result.sum = 0.0f;
		Result.Lnorm = (Params.Lnorm == 0 ? Path.K : Params.Lnorm);
		Result.anchor_path_indexes.clear();

		if (Path.K == 0 || Result.Lnorm == 0)
			return;

		ChooseAnchors(Path, Params.anchor_count, Params.min_anchor_sep,
			Result.anchor_path_indexes);

		Result.per_pair_score.assign(Path.K, 0.0f);
		Result.delta_Ang.assign(Path.K, 0.0f);
		Result.chosen_target.resize(Path.K);

		const float clip2 = Params.clip_Ang*Params.clip_Ang;
		const float d02 = Params.d0_Ang*Params.d0_Ang;

		float Sum = 0.0f;

		for (uint t = 0; t < Path.K; ++t)
			{
			if (Params.exclude_anchors_from_score &&
				IsAnchor(t, Result.anchor_path_indexes))
				continue;

			const uint i = Path.A[t];
			const uint j = Path.B[t];

			std::vector<seg_target_t> Targets;
			EnumTargets(j, DM_Y.m_L, Params, Targets);

			float best_delta2 = clip2;
			float best_score = 1.0f/(1.0f + clip2/d02);
			seg_target_t best_target;
			bool got_any = false;

			for (uint z = 0; z < Targets.size(); ++z)
				{
				float delta2;
				if (!EvalTargetDelta2(DM_X, DM_Y, Path, Result.anchor_path_indexes,
					i, Targets[z], t, clip2, delta2))
					continue;

				const float s = 1.0f/(1.0f + delta2/d02);
				if (!got_any || s > best_score)
					{
					got_any = true;
					best_score = s;
					best_delta2 = delta2;
					best_target = Targets[z];
					}
				}

			if (!got_any)
				{
				best_delta2 = clip2;
				best_score = 1.0f/(1.0f + clip2/d02);
				best_target.kind = SEG_TARGET_RESIDUE;
				best_target.j0 = (uint16_t) j;
				best_target.j1 = (uint16_t) j;
				best_target.lambda_num = 0;
				best_target.lambda_den = 1;
				}

			Result.per_pair_score[t] = best_score;
			Result.delta_Ang[t] = sqrtf(best_delta2);
			Result.chosen_target[t] = best_target;
			Sum += best_score;
			}

		Result.sum = Sum;
		Result.score = Sum/(float) Result.Lnorm;
	}

private:
	static bool IsAnchor(uint t, const std::vector<uint16_t>& Anchors)
	{
		for (uint i = 0; i < Anchors.size(); ++i)
			if (Anchors[i] == t)
				return true;
		return false;
	}

	static void ChooseAnchors(
		const align_path_t& Path,
		uint AnchorCount,
		uint MinSep,
		std::vector<uint16_t>& Anchors)
	{
		Anchors.clear();

		if (Path.K == 0 || AnchorCount == 0)
			return;

		if (Path.K <= AnchorCount)
			{
			for (uint t = 0; t < Path.K; ++t)
				Anchors.push_back((uint16_t) t);
			return;
			}

		uint first = 0;
		uint last = Path.K - 1;

		for (uint a = 0; a < AnchorCount; ++a)
			{
			uint t;
			if (AnchorCount == 1)
				t = (first + last)/2;
			else
				t = first + (uint)(((uint64_t) a*(last - first) + (AnchorCount - 1)/2)/
					(AnchorCount - 1));

			if (!Anchors.empty())
				{
				uint prev = Anchors.back();
				if (t <= prev)
					t = prev + 1;
				if (MinSep > 1 && t < prev + MinSep)
					t = prev + MinSep;
				}
			if (t > last)
				break;
			Anchors.push_back((uint16_t) t);
			}

		if (Anchors.empty())
			Anchors.push_back((uint16_t) ((first + last)/2));
	}

	static void AddTarget(std::vector<seg_target_t>& Targets,
		uint kind, uint j0, uint j1, uint lambda_num, uint lambda_den)
	{
		seg_target_t T;
		T.kind = (uint16_t) kind;
		T.j0 = (uint16_t) j0;
		T.j1 = (uint16_t) j1;
		T.lambda_num = (uint8_t) lambda_num;
		T.lambda_den = (uint8_t) lambda_den;
		Targets.push_back(T);
	}

	static void EnumTargets(
		uint j,
		uint LY,
		const rotfree_tm_segment_params_t& Params,
		std::vector<seg_target_t>& Targets)
	{
		Targets.clear();

		if (Params.allow_residue && j < LY)
			AddTarget(Targets, SEG_TARGET_RESIDUE, j, j, 0, 1);

		if (Params.allow_left_midpoint && j >= 1)
			AddTarget(Targets, SEG_TARGET_SEGMENT, j-1, j, 1, 2);

		if (Params.allow_right_midpoint && j + 1 < LY)
			AddTarget(Targets, SEG_TARGET_SEGMENT, j, j+1, 1, 2);

		if (Params.extra_lambda_den > 0 && !Params.extra_lambda_nums.empty())
			{
			for (uint k = 0; k < Params.extra_lambda_nums.size(); ++k)
				{
				uint num = Params.extra_lambda_nums[k];
				uint den = Params.extra_lambda_den;
				if (num == 0 || num >= den)
					continue;

				if (j >= 1)
					AddTarget(Targets, SEG_TARGET_SEGMENT, j-1, j, num, den);
				if (j + 1 < LY)
					AddTarget(Targets, SEG_TARGET_SEGMENT, j, j+1, num, den);
				}
			}
	}

	// -----------------------------------------------------------
	// Evaluate delta^2 for residue i in X against one candidate
	// target on Y backbone, using anchor-distance profiles.
	//
	// delta^2 = average_r min((dX - dYtarget)^2, clip2)
	//
	// Returns false if no anchors contributed.
	// -----------------------------------------------------------
	static bool EvalTargetDelta2(
		const flat_sidmx_t& DM_X,
		const flat_sidmx_t& DM_Y,
		const align_path_t& Path,
		const std::vector<uint16_t>& AnchorPathIndexes,
		uint i,
		const seg_target_t& Target,
		uint self_path_index,
		float clip2,
		float& delta2)
	{
		float sum_rho = 0.0f;
		uint n = 0;

		for (uint q = 0; q < AnchorPathIndexes.size(); ++q)
			{
			const uint a = AnchorPathIndexes[q];
			if (a == self_path_index)
				continue;

			const uint ia = Path.A[a];
			const uint ja = Path.B[a];

			if (!DM_X.HasPair(i, ia))
				continue;

			float dX = sid_to_ang(DM_X.Get(i, ia));
			float dY;
			if (!GetTargetDistAng(DM_Y, Target, ja, dY))
				continue;

			float diff = dX - dY;
			float e2 = diff*diff;
			if (e2 > clip2)
				e2 = clip2;

			sum_rho += e2;
			++n;
			}

		if (n == 0)
			return false;

		delta2 = sum_rho/(float) n;
		return true;
	}

	// -----------------------------------------------------------
	// Distance in Angstroms from anchor residue k to target point.
	//
	// Residue target:
	//   d = d(j,k)
	//
	// Segment target:
	//   d^2 = (1-lam)*d^2(j0,k) + lam*d^2(j1,k) - lam*(1-lam)*d^2(j0,j1)
	//
	// Requires all needed sid values to exist in the banded matrix.
	// -----------------------------------------------------------
	static bool GetTargetDistAng(
		const flat_sidmx_t& DM_Y,
		const seg_target_t& Target,
		uint k,
		float& dAng)
	{
		if (Target.kind == SEG_TARGET_RESIDUE)
			{
			uint j = Target.j0;
			if (j == k)
				{
				dAng = 0.0f;
				return true;
				}
			if (!DM_Y.HasPair(j, k))
				return false;

			dAng = sid_to_ang(DM_Y.Get(j, k));
			return true;
			}

		const uint j0 = Target.j0;
		const uint j1 = Target.j1;
		const float lam = (float) Target.lambda_num/(float) Target.lambda_den;

		float d2_0k = 0.0f;
		float d2_1k = 0.0f;
		float d2_01 = 0.0f;

		if (j0 == k)
			d2_0k = 0.0f;
		else
			{
			if (!DM_Y.HasPair(j0, k))
				return false;
			d2_0k = sid_to_ang2(DM_Y.Get(j0, k));
			}

		if (j1 == k)
			d2_1k = 0.0f;
		else
			{
			if (!DM_Y.HasPair(j1, k))
				return false;
			d2_1k = sid_to_ang2(DM_Y.Get(j1, k));
			}

		if (j0 == j1)
			d2_01 = 0.0f;
		else
			{
			if (!DM_Y.HasPair(j0, j1))
				return false;
			d2_01 = sid_to_ang2(DM_Y.Get(j0, j1));
			}

		float d2 =
			(1.0f - lam)*d2_0k +
			lam*d2_1k -
			lam*(1.0f - lam)*d2_01;

		if (d2 < 0.0f)
			d2 = 0.0f;

		dAng = sqrtf(d2);
		return true;
	}
};
