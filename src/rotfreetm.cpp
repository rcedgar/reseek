#include "myutils.h"
#include "rotfreetm.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "flat_bench_struct_feature.h"

void cmd_rotfreetm()
	{
//	const sid_t *DistA, uint LA,
//	const sid_t *DistB, uint LB,
//	uint M,
//	const uint16_t *AliA,
//	const uint16_t *AliB,
//	uint K)
//{
//	flat_sidmx_t DM_A(DistA, LA, M);
//	flat_sidmx_t DM_B(DistB, LB, M);
//
//	align_path_t Path;
//	Path.A = AliA;
//	Path.B = AliB;
//	Path.K = K;
//
//	rotfree_tm_params_t P;
//	P.anchor_count = 16;
//	P.min_anchor_sep = 8;
//	P.clip_Ang = 8.0f;
//	P.d0_Ang = 3.0f;
//	P.Lnorm = K;
//
//	rotfree_tm_scorer Scorer;
//	float Score = Scorer.ScoreAlignment(DM_A, DM_B, Path, P);
	}
	
float flat_bench_struct_feature::get_rotfreetm(uint idxQ, uint idxT,
	const flat_aligner &fa) const
	{
	asserta(idxQ < m_distmxs.size());
	asserta(idxT < m_distmxs.size());

	const sid_t *distmxQ = m_distmxs[idxQ];
	const sid_t *distmxT = m_distmxs[idxT];

	uint loQ = fa.m_loQ;
	uint loT = fa.m_loT;
	uint LQ = fa.m_LQ;
	uint LT = fa.m_LT;

	string path;
	uint nmatch = fa.get_path_str(path);
	uint ncol = uint(path.size());
	const uint M = flat_params::m_distmx_bandwidth;
	flat_sidmx_t DM_Q(distmxQ, LQ, M);
	flat_sidmx_t DM_T(distmxT, LT, M);

	vector<uint32_t> posQs;
	vector<uint32_t> posTs;
	path2posvecs(path, loQ, LQ, loT, LT, posQs, posTs);

	align_path_t Path;
	Path.A = posQs.data();
	Path.B = posTs.data();
	Path.K = K;

	rotfree_tm_params_t P;
	P.anchor_count = 16;
	P.min_anchor_sep = 8;
	P.clip_Ang = 8.0f;
	P.d0_Ang = 3.0f;
	P.Lnorm = K;

	rotfree_tm_scorer Scorer;
	float score = Scorer.ScoreAlignment(DM_Q, DM_T, Path, P);
	return score;
	}
