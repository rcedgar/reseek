#include "myutils.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "flat_features.h"
#include "flat_alignx.h"

// float L = float(LA + LB)/2;
//if (m_SelfRevScoreA != FLT_MAX && m_SelfRevScoreB != FLT_MAX)
//	RevDPScore = (m_SelfRevScoreA + m_SelfRevScoreB)/2;
//m_NewTestStatisticA = DSSParams::m_lddtw*LDDT;
//m_NewTestStatisticA += (DSSParams::m_dpw*m_AlnFwdScore -
//	DSSParams::m_revtsw*RevDPScore)/(L + DSSParams::m_ladd);
#if 1
static float oldts(
	const flat_aligner &fa,
	const uint8_t *profQ,
	const uint8_t *profT,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	float selfT,
	float selfQ,
	uint M)
	{
	asserta(distmxQ != 0);
	asserta(distmxT != 0);
	float LQ = (float) fa.m_LQ;
	float LT = (float) fa.m_LT;
	asserta(selfT != FLT_MAX && selfQ != FLT_MAX);
	float RevDPScore = (selfT + selfQ)/2;
	float LDDT = flat_getlddt_old(fa, distmxQ, distmxT, M);
	float AlnFwdScore = fa.m_score;
	float L = (LQ + LT)/2;
	float TS = flat_params::m_oldts_lddtw*LDDT;
	TS += (flat_params::m_oldts_dpw*AlnFwdScore - 
			flat_params::m_oldts_revtsw*RevDPScore)/(L + flat_params::m_oldts_ladd);
	return TS;
	}
#else
static float oldts(
	const flat_aligner &fa,
	const uint8_t *profQ,
	const uint8_t *profT,
	const sid_t *distmxQ,
	const sid_t *distmxT,
	float selfT,
	float selfQ,
	uint M)
	{
	asserta(distmxQ != 0);
	asserta(distmxT != 0);
	float LQ = (float) fa.m_LQ;
	float LT = (float) fa.m_LT;
	asserta(selfT != FLT_MAX && selfQ != FLT_MAX);
	float RevDPScore = (selfT + selfQ)/2;
	float LDDT = flat_getlddt_old(fa, distmxQ, distmxT, M);
	float AlnFwdScore = fa.m_score;
	float L = (LQ + LT)/2;
	uint match_count = fa.get_match_count();
	float TS = flat_params::m_oldts_lddtw*LDDT;
	TS += (flat_params::m_oldts_dpw*AlnFwdScore - 
			flat_params::m_oldts_revtsw*RevDPScore)/match_count;
	return TS;
	}
#endif

float flat_alignx::alignx(
	const flat_aligner &fa,
	const uint8_t *profQ, const uint8_t *profT,
	const sid_t *distmxQ, const sid_t *distmxT, uint M,
	float selfT, float selfQ)
	{
	float Score = fa.m_score;
	if (flat_params::m_oldts)
		Score = oldts(fa, profQ, profT, distmxQ, distmxT, selfT, selfQ, M);

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
