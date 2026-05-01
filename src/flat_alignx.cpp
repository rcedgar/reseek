#include "myutils.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "flat_features.h"
#include "flat_alignx.h"

float flat_alignx::alignx(
	const flat_aligner &fa,
	const uint8_t *profQ, const uint8_t *profT,
	const sid_t *distmxQ, const sid_t *distmxT,
	float selfT, float selfQ)
	{
	float Score = fa.m_score;

	if (flat_params::m_self_w > 0)
		{
		asserta(selfT != FLT_MAX && selfQ != FLT_MAX);
		Score -= flat_params::m_self_w*(selfT + selfQ)/2;
		}

	if (flat_params::m_lddt_w > 0)
		{
		asserta(distmxQ != 0 && distmxT != 0);
		float lddt = flat_getlddt_muscle_some_floats4(
			fa, distmxQ, distmxT);
		Score += flat_params::m_lddt_w*lddt*500;
		}

	if (flat_params::m_lddtx_w > 0)
		{
		float LQ = (float) fa.m_LQ;
		float LT = (float) fa.m_LT;
		float L = (LQ + LT)/2.0f + 50;

		uint ncol = uint(fa.m_ncol);
		float Lfactor = float(ncol)/L;

		asserta(distmxQ != 0 && distmxT != 0);
		float lddt = flat_getlddt_muscle_some_floats4(
			fa, distmxQ, distmxT);
		Score += flat_params::m_lddtx_w*lddt*500*Lfactor;
		}

	if (flat_params::m_lddtpow_w > 0)
		{
		float LQ = (float) fa.m_LQ;
		float LT = (float) fa.m_LT;
		float L = (LQ + LT)/2.0f + 50;

		uint ncol = uint(fa.m_ncol);
		float Lfactor = float(ncol)/L;

		asserta(distmxQ != 0 && distmxT != 0);
		float lddt = flat_getlddt_muscle_some_floats4(
			fa, distmxQ, distmxT);
		float maxL = max(LT, LQ) - 20.0f;
		if (maxL < 80)
			maxL = 80;
		string path;
		uint nmatch = fa.get_path_str(path);
		float lddtpow = lddt*nmatch*2.0f/powf(maxL, 0.5);
		Score += flat_params::m_lddtpow_w*lddtpow;
		}

	if (flat_params::m_dali_w > 0)
		{
		asserta(distmxQ != 0 && distmxT != 0);
		float dali = flat_get_dali3(fa, distmxQ, distmxT);
		Score += flat_params::m_dali_w*dali*10;
		}

	if (flat_params::m_dalix_w > 0)
		{
		asserta(distmxQ != 0 && distmxT != 0);
		uint ncol = uint(fa.m_ncol);
		float *colscores = myalloc(float, ncol);
		float dalix = flat_get_dalix3(fa, distmxQ, distmxT, colscores);
		Score += flat_params::m_dalix_w*dalix*10;
		myfree(colscores);
		}

	if (flat_params::m_entropy_w > 0)
		{
		asserta(profQ != 0 && profT != 0);
		float entropy = flat_get_entropy2(
			fa, profQ, profT,
			flat_features::m_nfeat,
			flat_features::m_entropyfi);
		Score += flat_params::m_entropy_w*entropy/250;
		}
	return Score;
	}
