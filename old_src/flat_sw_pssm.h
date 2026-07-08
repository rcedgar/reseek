#pragma once

template<uint32_t NFEAT>
float smith_waterman_affine_flat_pssm_fixed_nfeat(
	const uint8_t * __restrict profA,
	uint32_t LA,
	uint32_t LB,
	const uint32_t * __restrict feature_block_offsets,
	const float * __restrict pssm,
	float gap_open,
	float gap_extend,
	float * __restrict H,
	float * __restrict F)
	{
	static_assert(NFEAT >= 1, "NFEAT must be >= 1");
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
		const float * __restrict rows[NFEAT];

		// Select one PSSM row per feature for this i.
		for (uint32_t fi = 0; fi < NFEAT; ++fi)
			{
			const uint8_t codeA_i = profA[size_t(fi)*LA + i];
			const float * __restrict pssm_fi =
				pssm + size_t(feature_block_offsets[fi])*LB;
			rows[fi] = pssm_fi + size_t(codeA_i)*LB;
			}

		float H_left = 0.0f;   // H(i,j-1)
		float H_diag = 0.0f;   // H(i-1,j-1)
		float E = NEG_INF;     // horizontal gap state

		for (uint32_t j = 0; j < LB; ++j)
			{
			const float H_up = H[j];

			float sub = 0.0f;
			for (uint32_t fi = 0; fi < NFEAT; ++fi)
				sub += rows[fi][j];

			const float e_open = H_left - gap_open;
			const float e_ext = E - gap_extend;
			E = (e_open > e_ext ? e_open : e_ext);

			const float f_open = H_up - gap_open;
			const float f_ext = F[j] - gap_extend;
			F[j] = (f_open > f_ext ? f_open : f_ext);

			float h = H_diag + sub;
			if (E > h)
				h = E;
			if (F[j] > h)
				h = F[j];
			if (h < 0.0f)
				h = 0.0f;

			H_diag = H_up;
			H[j] = h;
			H_left = h;

			if (h > best_score)
				best_score = h;
			}
		}

	return best_score;
	}
