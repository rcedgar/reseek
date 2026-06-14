#include "myutils.h"
#include "reseeker.h"
#include "flat_nu_aligner.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "hitdata.h"

static double s_min_ts_fold = 30;//TODO param

struct mega_ts_result
	{
	float mega_fwd_score = 0;
	float mega_rev_score = 0;
	float mega_self_score = 0;
	float lddt = 0;
	float dali = 0;
	float TS = 0;
	};

static float get_target_mega_self_rev_score(
	float &cached_score,
	const flat_params &params,
	const uint8_t *target_mega_prof,
	uint LT,
	const float *target_mega_pssm_rev,
	float *target_pssm_rev_buf,
	float *scratch_rows,
	const float **scratch_pssms)
	{
	if (cached_score != FLT_MAX)
		return cached_score;

	const float *pssm_rev = target_mega_pssm_rev;
	if (pssm_rev == 0)
		{
		fill_flat_pssm_reversed(
			target_mega_prof, LT,
			params.m_nfeat, params.m_alpha_sizes,
			params.m_feature_block_offsets,
			params.m_weighted_logoddsvec,
			target_pssm_rev_buf);
		pssm_rev = target_pssm_rev_buf;
		}

	cached_score = sw_flat_pssm_scoreonly(
		scratch_rows, scratch_pssms,
		target_mega_prof, LT, pssm_rev, LT,
		params.m_feature_block_offsets,
		params.m_nfeat,
		-params.m_open,
		-params.m_ext);
	return cached_score;
	}

static mega_ts_result compute_mega_ts(
	const flat_params &params,
	const uint8_t *target_mega_prof,
	uint LT,
	const float *query_mega_pssm,
	const float *query_mega_pssm_rev,
	uint LQ,
	float query_mega_self_rev_score,
	float &target_mega_self_rev_score,
	const float *target_mega_pssm_rev,
	float *target_pssm_rev_buf,
	int nu_rev_score,
	const sid_t *distmx_i,
	const sid_t *distmx_j,
	float *scratch_rows,
	uint8_t *TB,
	const float **scratch_pssms,
	char *path_buffer,
	uint *pos_is,
	uint *pos_js,
	uint *considered_vec,
	uint *preserved_vec,
	uint maxL)
	{
	mega_ts_result result;

	uint lo_i, lo_j, ncol;
	result.mega_fwd_score = sw_flat_pssm(
		scratch_rows, TB, scratch_pssms,
		target_mega_prof, LT, query_mega_pssm, LQ,
		params.m_feature_block_offsets,
		params.m_nfeat,
		-params.m_open,
		-params.m_ext,
		lo_i, lo_j, path_buffer, ncol);

	float target_self_rev = get_target_mega_self_rev_score(
		target_mega_self_rev_score,
		params, target_mega_prof, LT,
		target_mega_pssm_rev, target_pssm_rev_buf,
		scratch_rows, scratch_pssms);

	result.mega_rev_score = sw_flat_pssm_scoreonly(
		scratch_rows, scratch_pssms,
		target_mega_prof, LT, query_mega_pssm_rev, LQ,
		params.m_feature_block_offsets,
		params.m_nfeat,
		-params.m_open,
		-params.m_ext);

	result.mega_self_score =
		(target_self_rev + query_mega_self_rev_score)/2;

	const uint L_i = LT;
	const uint L_j = LQ;
	uint nmatch = path2posvecs3(path_buffer, ncol,
		lo_i, L_i, lo_j, L_j, pos_is, pos_js, maxL);

	result.lddt = flat_getlddt_muscle_some_floats(
		pos_is, L_i,
		pos_js, L_j,
		nmatch, distmx_i, distmx_j,
		considered_vec, preserved_vec);

	result.dali = flat_get_dali4(
		pos_is, L_i,
		pos_js, L_j,
		nmatch, distmx_i, distmx_j);

	asserta(params.m_lddtx_w == 0);
	asserta(params.m_dalix_w == 0);

	result.TS = 0;
	result.TS += result.mega_fwd_score;
	result.TS -= params.m_rev_w*result.mega_rev_score;
	result.TS -= params.m_self_w*result.mega_self_score;
	result.TS += params.m_nurev_w*nu_rev_score;
	result.TS += params.m_lddt_w*result.lddt*500;
	result.TS += params.m_dali_w*result.dali*10;

	return result;
	}

void reseeker::static_thread_body(uint threadidx)
	{
	const BCAData &dbbca = *m_dbbca;

	const int nu_open = flat_nu_aligner::m_open;
	const int nu_ext = flat_nu_aligner::m_ext;

	const flat_params &params_fold = *m_params_fold;
	const float revw = params_fold.m_nu_filter_rev_w;
	const float selfw = params_fold.m_nu_filter_self_w;

	uint workspace_bytes =
		parasail_nomalloc_sw_striped_profile_avx2_256_16_workspace_bytes(m_maxL);
	uint scratch_buffer_bytes = 2*m_maxL;

	uint8_t *workspace = myalloc(uint8_t, workspace_bytes);
	float *scratch_rows = myalloc(float, 2*m_maxL + 2);
	const float **scratch_pssms =
		myalloc(const float *, params_fold.m_nfeat);
	uint8_t *TB = myalloc(uint8_t, m_maxL*m_maxL);
	char *path_buffer = myalloc(char, 2*m_maxL);
	uint *pos_is = myalloc(uint, m_maxL);
	uint *pos_js = myalloc(uint, m_maxL);
	uint *considered_vec = myalloc(uint, m_maxL);
	uint *preserved_vec = myalloc(uint, m_maxL);
	uint8_t *scratch_buffer = myalloc(uint8_t, scratch_buffer_bytes);
	float *target_pssm_rev_buf = myalloc(
		float, m_maxL*params_fold.m_sum_alpha_sizes);
	chaq_vecs2 cv;
	chaq::alloc_chaq_vecs2(cv, m_maxL);
	hitdata hit;

	for (;;)
		{
		uint k = m_next++;
		if (k + 1 == m_ndbidxs || (k%100 == 0 && k > 0 && k + 1 < m_ndbidxs))
			{
			static mutex progress_lock;
			progress_lock.lock();
			ProgressStep(k, m_ndbidxs, "reseek");
			progress_lock.unlock();
			}
		if (k >= m_ndbidxs) return;

		const vector<uint> *ptr_qidxs = 0;
		uint dbidx = UINT_MAX;

		switch (m_mode)
			{
		case NF_all_vs_all:
			{
			dbidx = k;
			ptr_qidxs = &m_qidxs_all;
			break;
			}

		case NF_kappa:
			{
			dbidx = (*m_dbidxs)[k];
			unordered_map<uint, vector<uint> >::const_iterator iter =
				m_dbidx_to_qidxs->find(dbidx);
			asserta(iter != m_dbidx_to_qidxs->end());
			ptr_qidxs = &iter->second;
			break;
			}

		default: asserta(false);
			}

		const vector<uint> &qidxs = *ptr_qidxs;
		const uint nq = uint(qidxs.size());
		asserta(nq > 0);
		struct_data *target_data = dbbca.get_struct_data(
			params_fold, dbidx,
			&cv, scratch_buffer, scratch_buffer_bytes);
		const string &target_label = dbbca.m_Labels[dbidx];
		uint8_t *target_codeseq_nu = target_data->m_codeseq_nu;
		uint8_t *target_codeseq_nu_rev = target_data->m_codeseq_nu_rev;
		uint8_t *target_mega_prof = target_data->m_mega_prof;
		const sid_t *target_distmx = target_data->m_distmx;
		const uint LT = target_data->m_chain->get_length();

		parasail_profile_t *target_para_prof = target_data->m_parasail_prof;
		int target_self_rev_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
			target_para_prof, (const char *) target_codeseq_nu_rev, LT, nu_open, nu_ext,
			workspace, workspace_bytes);

		float target_mega_self_rev_score_fold = FLT_MAX;
		float target_mega_self_rev_score_sf = FLT_MAX;
		float target_mega_self_rev_score_fam = FLT_MAX;
		for (uint j = 0; j < nq; ++j)
			{
			++m_npair;
			uint qidx = qidxs[j];
			asserta(qidx < m_query_nchain);
			const string &query_label = (*m_ptr_query_labels)[qidx];
			parasail_profile_t *query_para_prof = m_query_parasail_profs[qidx];

			/////////////////////////////////////////////
			// Nu filter -- forward score
			/////////////////////////////////////////////
			int nu_fwd_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
				query_para_prof, (const char *) target_codeseq_nu, LT, nu_open, nu_ext,
				workspace, workspace_bytes);
			if (nu_fwd_score < params_fold.m_nu_filter_min_fwd_score)
				{
				++m_reject_fwd;
				continue;
				}

			parasail_profile_t *query_para_prof_rev = m_query_parasail_prof_revs[qidx];
			int nu_rev_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
				query_para_prof_rev, (const char *) target_codeseq_nu, LT, nu_open, nu_ext,
				workspace, workspace_bytes);

			float self_score = (target_self_rev_score + 
				m_query_self_rev_scores[qidx])/2.0f;

			/////////////////////////////////////////////
			// Nu filter -- combined score
			/////////////////////////////////////////////
			float combined_score = float(nu_fwd_score) - selfw*self_score - revw*nu_rev_score;
			if (combined_score < params_fold.m_nu_filter_min_combined_score)
				{
				++m_reject_cmb;
				continue;
				}
			++m_npass;

			const uint LQ = m_query_lengths[qidx];
			const sid_t *distmx_j = m_query_distmxs[qidx];

			mega_ts_result fold_ts = compute_mega_ts(
				params_fold,
				target_mega_prof, LT,
				m_query_mega_pssms_fold[qidx],
				m_query_mega_pssm_revs_fold[qidx],
				LQ,
				m_query_mega_self_rev_scores_fold[qidx],
				target_mega_self_rev_score_fold,
				target_data->m_mega_pssm_rev,
				target_pssm_rev_buf,
				nu_rev_score,
				target_distmx, distmx_j,
				scratch_rows, TB, scratch_pssms,
				path_buffer, pos_is, pos_js,
				considered_vec, preserved_vec,
				m_maxL);

			if (fold_ts.mega_fwd_score < params_fold.m_mega_filter_min_fwd)
				{
				++m_reject_mega_fwd;
				continue;
				}

			if (fold_ts.TS < s_min_ts_fold)
				{
				++m_reject_min_ts;
				continue;
				}
			++m_accept_min_ts;

			float TS_sf = 0;
			if (m_params_sf != 0)
				{
				mega_ts_result sf_ts = compute_mega_ts(
					*m_params_sf,
					target_mega_prof, LT,
					m_query_mega_pssms_sf[qidx],
					m_query_mega_pssm_revs_sf[qidx],
					LQ,
					m_query_mega_self_rev_scores_sf[qidx],
					target_mega_self_rev_score_sf,
					0,
					target_pssm_rev_buf,
					nu_rev_score,
					target_distmx, distmx_j,
					scratch_rows, TB, scratch_pssms,
					path_buffer, pos_is, pos_js,
					considered_vec, preserved_vec,
					m_maxL);
				TS_sf = sf_ts.TS;
				}

			float TS_fam = 0;
			if (m_params_fam != 0)
				{
				mega_ts_result fam_ts = compute_mega_ts(
					*m_params_fam,
					target_mega_prof, LT,
					m_query_mega_pssms_fam[qidx],
					m_query_mega_pssm_revs_fam[qidx],
					LQ,
					m_query_mega_self_rev_scores_fam[qidx],
					target_mega_self_rev_score_fam,
					0,
					target_pssm_rev_buf,
					nu_rev_score,
					target_distmx, distmx_j,
					scratch_rows, TB, scratch_pssms,
					path_buffer, pos_is, pos_js,
					considered_vec, preserved_vec,
					m_maxL);
				TS_fam = fam_ts.TS;
				}

			hit.reset();
			hit.query = m_ptr_query_chains[qidx];
			hit.target = target_data->m_chain;
			hit.nu_fwd_score = float(nu_fwd_score);
			hit.nu_rev_score = float(nu_rev_score);
			hit.mega_fwd_score = fold_ts.mega_fwd_score;
			hit.mega_rev_score = fold_ts.mega_rev_score;
			hit.lddt = fold_ts.lddt;
			hit.dali = fold_ts.dali;
			hit.TS_fold = fold_ts.TS;
			hit.TS_sf = TS_sf;
			hit.TS_fam = TS_fam;

			asserta(m_fhit);
			if (m_fhit)
				{
				string str;
				str = query_label;
				str += "\t" + target_label;
				Psa(str, "\t%.3g", hit.TS_fold);
				if (m_params_sf != 0)
					Psa(str, "\t%.3g", hit.TS_sf);
				if (m_params_fam != 0)
					Psa(str, "\t%.3g", hit.TS_fam);
				//Psa(str, "\t%.3g", float(nu_fwd_score));
				//Psa(str, "\t%.3g", float(nu_rev_score));
				//Psa(str, "\t%.3g", fold_ts.mega_fwd_score);
				//Psa(str, "\t%.3g", fold_ts.mega_rev_score);
				//Psa(str, "\t%.3g", fold_ts.lddt);
				//Psa(str, "\t%.3g", fold_ts.dali);
				str += "\n";
				fputs(str.c_str(), m_fhit);
				}
			}

		struct_data::free_struct_data(target_data);
		}
	}
