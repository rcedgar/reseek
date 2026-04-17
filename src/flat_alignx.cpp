#include "myutils.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "flat_features.h"
#include "flat_alignx.h"

float flat_alignx::alignx(
	const flat_aligner &fa,
	const uint8_t *profQ,
	const uint8_t *profT,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	float selfT,
	float selfQ,
	uint M)
	{
	float Score = 0;
	if (flat_params::m_rev_w > 0)
		{
		asserta(fa.m_reverse_score_set);
		Score -= flat_params::m_rev_w*fa.m_reverse_score;
		}

	if (flat_params::m_self_w > 0)
		{
		asserta(selfT != FLT_MAX && selfQ != FLT_MAX);
		Score -= flat_params::m_self_w*(selfT + selfQ);
		}

	if (flat_params::m_lddt_w > 0)
		{
		asserta(distmxQ != 0 && distmxT != 0);
		float lddt = flat_getlddt_muscle_some_floats4(
			fa, distmxQ, distmxT, M);
		Score += flat_params::m_lddt_w*lddt/4;
		}

	if (flat_params::m_dali_w > 0)
		{
		asserta(distmxQ != 0 && distmxT != 0);
		uint n = uint(fa.m_ncol);
		float dali = 
			flat_get_dali3(fa, distmxQ, distmxT, M);
		Score += flat_params::m_dali_w*dali*10;
		}

	if (flat_params::m_dalix_w > 0)
		{
		asserta(distmxQ != 0 && distmxT != 0);
		uint n = uint(fa.m_ncol);
		float *colscores = myalloc(float, n);
		float dalix = 
			flat_get_dalix3(fa, distmxQ, distmxT, M, colscores);
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
