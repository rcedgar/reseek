#include "myutils.h"
#include "reseeker.h"
#include "reseeker_thread_body_impl.h"
#include "flat_nu_aligner.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "hitdata.h"
#include "reseek_hit_sink.h"
#include <algorithm>

static void resolve_target(uint k,
	uint &dbidx,
	const vector<uint> **ptr_qidxs,
	const vector<uint> **ptr_diagscores)
	{
	*ptr_qidxs = 0;
	*ptr_diagscores = 0;
	dbidx = UINT_MAX;

	switch (reseeker::m_mode)
		{
	case NF_all_vs_all:
		dbidx = k;
		*ptr_qidxs = &reseeker::m_qidxs_all;
		break;

	case NF_kappa:
		dbidx = (*reseeker::m_dbidxs)[k];
		{
		unordered_map<uint, vector<uint> >::const_iterator iter =
			reseeker::m_dbidx_to_qidxs->find(dbidx);
		asserta(iter != reseeker::m_dbidx_to_qidxs->end());
		*ptr_qidxs = &iter->second;
		unordered_map<uint, vector<uint> >::const_iterator iter_diag =
			reseeker::m_dbidx_to_diagscores->find(dbidx);
		asserta(iter_diag != reseeker::m_dbidx_to_diagscores->end());
		*ptr_diagscores = &iter_diag->second;
		asserta((*ptr_diagscores)->size() == (*ptr_qidxs)->size());
		}
		break;

	default:
		asserta(false);
		}
	}

static bool nu_pass_gt(const nu_pass &a, const nu_pass &b)
	{
	return a.nu_combined_score > b.nu_combined_score;
	}

static void finalize_nu_passes(vector<nu_pass> &passes, bool use_nusort)
	{
	if (!use_nusort)
		return;
	const uint max_nu_accepts = flat_params::m_max_nu_filter_accepts;
	const uint n = uint(passes.size());
	if (n <= max_nu_accepts)
		return;
	std::nth_element(passes.begin(), passes.begin() + max_nu_accepts, passes.end(),
		nu_pass_gt);
	std::sort(passes.begin(), passes.begin() + max_nu_accepts, nu_pass_gt);
	passes.resize(max_nu_accepts);
	}

static bool score_mega_hit(
	const flat_params &params,
	float *scratch_rows,
	const float **scratch_pssms,
	uint8_t *TB,
	char *path_buffer,
	uint *pos_is,
	uint *pos_js,
	uint *considered_vec,
	uint *preserved_vec,
	struct_data *target_data,
	const chain_slice &target_slice,
	float &target_mega_self_rev_score,
	uint qidx,
	int nu_fwd_score,
	int nu_rev_score,
	float nu_combined_score,
	uint kappa_diag_score,
	const string &query_label,
	const string &target_label,
	reseek_hit &hit)
	{
	const uint LQ = reseeker::m_query_lengths[qidx];
	const float *query_mega_pssm = reseeker::m_query_mega_pssms[qidx];
	const float *query_mega_pssm_rev = reseeker::m_query_mega_pssm_revs[qidx];
	uint8_t *target_mega_prof = target_data->m_mega_prof;
	const sid_t *target_distmx = target_data->m_distmx;
	const uint LT = target_data->m_chain->get_length();

	uint lo_i, lo_j, ncol;
	float mega_fwd_score = sw_flat_pssm(
		scratch_rows, TB, scratch_pssms,
		target_mega_prof, LT, query_mega_pssm, LQ,
		params.m_feature_block_offsets,
		params.m_nfeat,
		-params.m_open,
		-params.m_ext,
		lo_i, lo_j, path_buffer, ncol);
	const uint L_i = LT;
	const uint L_j = LQ;
	const sid_t *distmx_i = target_distmx;
	const sid_t *distmx_j = reseeker::m_query_distmxs[qidx];

	if (mega_fwd_score < params.m_mega_filter_min_fwd)
		{
		++reseeker::m_reject_mega_fwd;
		return false;
		}

	if (target_mega_self_rev_score == FLT_MAX)
		{
		const float *target_mega_pssm_rev = target_data->m_mega_pssm_rev;
		target_mega_self_rev_score = sw_flat_pssm_scoreonly(
			scratch_rows, scratch_pssms,
			target_mega_prof, LT, target_mega_pssm_rev, LT,
			params.m_feature_block_offsets,
			params.m_nfeat,
			-params.m_open,
			-params.m_ext);
		}

	float mega_rev_score = sw_flat_pssm_scoreonly(
		scratch_rows, scratch_pssms,
		target_mega_prof, LT, query_mega_pssm_rev, LQ,
		params.m_feature_block_offsets,
		params.m_nfeat,
		-params.m_open,
		-params.m_ext);

	float query_mega_self_rev_score =
		reseeker::m_query_mega_self_rev_scores[qidx];

	uint nmatch = path2posvecs3(path_buffer, ncol,
		lo_i, L_i, lo_j, L_j, pos_is, pos_js, flat_params::m_maxL);

	float lddt = flat_getlddt_muscle_some_floats(
		pos_is, L_i,
		pos_js, L_j,
		nmatch, distmx_i, distmx_j,
		considered_vec, preserved_vec);

	float dali = flat_get_dali4(
		pos_is, L_i,
		pos_js, L_j,
		nmatch, distmx_i, distmx_j);

	asserta(params.m_lddtx_w == 0);
	asserta(params.m_dalix_w == 0);

	asserta(target_mega_self_rev_score != FLT_MAX);
	float mega_self_score =
		(target_mega_self_rev_score + query_mega_self_rev_score)/2;

	float TS = 0;
	TS += mega_fwd_score;
	TS -= params.m_rev_w*mega_rev_score;
	TS -= params.m_self_w*mega_self_score;
	TS += params.m_nurev_w*nu_rev_score;
	TS += params.m_lddt_w*lddt*500;
	TS += params.m_dali_w*dali*10;

	if (TS < reseeker::m_mints)
		{
		++reseeker::m_reject_min_ts;
		return false;
		}
	++reseeker::m_accept_min_ts;

	hit.reset();
	hit.query = reseeker::m_query_slices[qidx];
	hit.target = target_slice;
	hit.parent_label_q = query_label;
	hit.parent_label_t = target_label;
	hit.set_fwd_path(path_buffer, ncol, lo_j, lo_i);
	reseek_hit_fill_aa_cigar(hit,
		reseeker::m_ptr_query_chains[qidx], target_data->m_chain);
	hit.nu_fwd_score = float(nu_fwd_score);
	hit.nu_rev_score = float(nu_rev_score);
	hit.nu_combined_score = nu_combined_score;
	hit.mega_fwd_score = mega_fwd_score;
	hit.mega_rev_score = mega_rev_score;
	hit.lddt = lddt;
	hit.dali = dali;
	hit.TS = TS;
	hit.kappa_diag_score = kappa_diag_score;
	return true;
	}

void reseeker_thread_body_impl(uint threadidx, bool use_nusort)
	{
	(void) threadidx;
	const BCAData &dbbca = *reseeker::m_dbbca;

	const int nu_open = flat_nu_aligner::m_open;
	const int nu_ext = flat_nu_aligner::m_ext;

	const flat_params &params = *reseeker::m_params;
	const float revw = reseeker::m_params->m_nu_filter_rev_w;
	const float selfw = reseeker::m_params->m_nu_filter_self_w;
	const bool nu_only = reseeker::m_params->m_nu_only;

	uint workspace_bytes =
		parasail_nomalloc_sw_striped_profile_avx2_256_16_workspace_bytes(flat_params::m_maxL);
	uint scratch_buffer_bytes = 2*flat_params::m_maxL;

	uint8_t *workspace = myalloc(uint8_t, workspace_bytes);
	float *scratch_rows = myalloc(float, 2*flat_params::m_maxL + 2);
	const float **scratch_pssms = myalloc(const float *, reseeker::m_params->m_nfeat);
	uint8_t *TB = myalloc(uint8_t, flat_params::m_maxL*flat_params::m_maxL);
	char *path_buffer = myalloc(char, 2*flat_params::m_maxL);
	uint *pos_is = myalloc(uint, flat_params::m_maxL);
	uint *pos_js = myalloc(uint, flat_params::m_maxL);
	uint *considered_vec = myalloc(uint, flat_params::m_maxL);
	uint *preserved_vec = myalloc(uint, flat_params::m_maxL);
	uint8_t *scratch_buffer = myalloc(uint8_t, scratch_buffer_bytes);
	chaq_vecs2 cv;
	chaq::alloc_chaq_vecs2(cv, flat_params::m_maxL);
	reseek_hit hit;

	for (;;)
		{
		uint k = reseeker::m_next++;
		if (k + 1 == reseeker::m_ndbidxs ||
			(k%100 == 0 && k > 0 && k + 1 < reseeker::m_ndbidxs))
			{
			static mutex progress_lock;
			progress_lock.lock();
			ProgressStep(k, reseeker::m_ndbidxs, "reseek");
			progress_lock.unlock();
			}
		if (k >= reseeker::m_ndbidxs)
			break;

		const vector<uint> *ptr_qidxs = 0;
		const vector<uint> *ptr_diagscores = 0;
		uint dbidx = UINT_MAX;
		resolve_target(k, dbidx, &ptr_qidxs, &ptr_diagscores);

		const vector<uint> &qidxs = *ptr_qidxs;
		const uint nq = uint(qidxs.size());
		asserta(nq > 0);
		struct_data *target_data = dbbca.get_struct_data(
			params, dbidx,
			&cv, scratch_buffer, scratch_buffer_bytes);
		const string &target_label = dbbca.m_Labels[dbidx];
		uint8_t *target_codeseq_nu = target_data->m_codeseq_nu;
		uint8_t *target_codeseq_nu_rev = target_data->m_codeseq_nu_rev;
		const uint LT = target_data->m_chain->get_length();
		chain_slice target_slice = chain_slice_identity(dbidx, LT);

		parasail_profile_t *target_para_prof = target_data->m_parasail_prof;
		int target_self_rev_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
			target_para_prof, (const char *) target_codeseq_nu_rev, LT,
			nu_open, nu_ext, workspace, workspace_bytes);

		vector<nu_pass> nu_passes;
		nu_passes.reserve(nq);

		for (uint j = 0; j < nq; ++j)
			{
			++reseeker::m_npair;
			uint qidx = qidxs[j];
			asserta(qidx < reseeker::m_query_nchain);
			parasail_profile_t *query_para_prof = reseeker::m_query_parasail_profs[qidx];

			int nu_fwd_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
				query_para_prof, (const char *) target_codeseq_nu, LT,
				nu_open, nu_ext, workspace, workspace_bytes);
			if (nu_fwd_score < params.m_nu_filter_min_fwd_score)
				{
				++reseeker::m_reject_fwd;
				continue;
				}

			parasail_profile_t *query_para_prof_rev =
				reseeker::m_query_parasail_prof_revs[qidx];
			int nu_rev_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
				query_para_prof_rev, (const char *) target_codeseq_nu, LT,
				nu_open, nu_ext, workspace, workspace_bytes);

			float self_score = (target_self_rev_score +
				reseeker::m_query_self_rev_scores[qidx])/2.0f;

			float nu_combined_score = float(nu_fwd_score) -
				selfw*self_score - revw*float(nu_rev_score);
			if (nu_combined_score < reseeker::m_params->m_nu_filter_min_combined_score)
				{
				++reseeker::m_nu_reject_cmb;
				continue;
				}
			++reseeker::m_npass;

			const string &query_label = (*reseeker::m_ptr_query_labels)[qidx];

			if (nu_only)
				{
				hit.reset();
				hit.query = reseeker::m_query_slices[qidx];
				hit.target = target_slice;
				hit.parent_label_q = query_label;
				hit.parent_label_t = target_label;
				hit.nu_combined_score = nu_combined_score;
				reseek_hit_sink_submit(hit, true);
				continue;
				}

			nu_pass p;
			p.qidx = qidx;
			p.j = j;
			p.nu_fwd_score = nu_fwd_score;
			p.nu_rev_score = nu_rev_score;
			p.nu_combined_score = nu_combined_score;
			nu_passes.push_back(p);
			}

		if (nu_only)
			{
			struct_data::free_struct_data(target_data);
			continue;
			}

		finalize_nu_passes(nu_passes, use_nusort);
		const uint n_top = uint(nu_passes.size());

		float target_mega_self_rev_score = FLT_MAX;
		for (uint i = 0; i < n_top; ++i)
			{
			const nu_pass &p = nu_passes[i];
			uint qidx = p.qidx;
			uint kappa_diag_score = 0;
			if (reseeker::m_mode == NF_kappa)
				kappa_diag_score = (*ptr_diagscores)[p.j];
			const string &query_label = (*reseeker::m_ptr_query_labels)[qidx];

			if (!score_mega_hit(params,
				scratch_rows, scratch_pssms, TB, path_buffer,
				pos_is, pos_js, considered_vec, preserved_vec,
				target_data, target_slice, target_mega_self_rev_score,
				qidx, p.nu_fwd_score, p.nu_rev_score, p.nu_combined_score,
				kappa_diag_score, query_label, target_label, hit))
				continue;

			asserta(reseeker::m_fhit);
			reseek_hit_sink_submit(hit, false);
			}

		struct_data::free_struct_data(target_data);
		}

	reseek_hit_sink_thread_end();

	chaq::free_chaq_vecs2(cv);
	myfree(scratch_buffer);
	myfree(preserved_vec);
	myfree(considered_vec);
	myfree(pos_js);
	myfree(pos_is);
	myfree(path_buffer);
	myfree(TB);
	myfree(scratch_pssms);
	myfree(scratch_rows);
	myfree(workspace);
	}

void reseeker::static_thread_body(uint threadidx)
	{
	reseeker_thread_body_impl(threadidx, false);
	}
