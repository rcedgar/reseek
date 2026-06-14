#include "myutils.h"
#include "reseeker.h"
#include "flat_nu_aligner.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "hitdata.h"

static double s_min_ts_fold = 30;//TODO param

void reseeker::static_thread_body(uint threadidx)
	{
	const BCAData &dbbca = *m_dbbca;

	const int nu_open = flat_nu_aligner::m_open;
	const int nu_ext = flat_nu_aligner::m_ext;

	const flat_params &params = *m_params;
	const float revw = m_params->m_nu_filter_rev_w;
	const float selfw = m_params->m_nu_filter_self_w;

	uint workspace_bytes =
		parasail_nomalloc_sw_striped_profile_avx2_256_16_workspace_bytes(m_maxL);
	uint scratch_buffer_bytes = 2*m_maxL;

	uint8_t *workspace = myalloc(uint8_t, workspace_bytes);
	float *scratch_rows = myalloc(float, 2*m_maxL + 2);
	const float **scratch_pssms = myalloc(const float *, m_params->m_nfeat);
	uint8_t *TB = myalloc(uint8_t, m_maxL*m_maxL);
	char *path_buffer = myalloc(char, 2*m_maxL);
	uint *pos_is = myalloc(uint, m_maxL);
	uint *pos_js = myalloc(uint, m_maxL);
	uint *considered_vec = myalloc(uint, m_maxL);
	uint *preserved_vec = myalloc(uint, m_maxL);
	uint8_t *scratch_buffer = myalloc(uint8_t, scratch_buffer_bytes);
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
			params, dbidx,
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

		float target_mega_self_rev_score = FLT_MAX; // calculate only if needed
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
			if (nu_fwd_score < params.m_nu_filter_min_fwd_score)
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
			if (combined_score < m_params->m_nu_filter_min_combined_score)
				{
				++m_reject_cmb;
				continue;
				}
			++m_npass;

			const uint LQ = m_query_lengths[qidx];
			const float *query_mega_pssm = m_query_mega_pssms[qidx];
			const float *query_mega_pssm_rev = m_query_mega_pssm_revs[qidx];

			/////////////////////////////////////////////
			// Mega forward score
			/////////////////////////////////////////////
			uint lo_i, lo_j, ncol;
			float mega_fwd_score = sw_flat_pssm(
				scratch_rows, TB, scratch_pssms,
				target_mega_prof, LT, query_mega_pssm, LQ,
				m_params->m_feature_block_offsets,
				m_params->m_nfeat,
				-m_params->m_open,
				-m_params->m_ext,
				lo_i, lo_j, path_buffer, ncol);
			const uint L_i = LT;
			const uint L_j = LQ;
			const sid_t *distmx_i = target_distmx;
			const sid_t *distmx_j = m_query_distmxs[qidx];

			/////////////////////////////////////////////
			// Reject if low mega forward score
			/////////////////////////////////////////////
			if (mega_fwd_score < m_params->m_mega_filter_min_fwd)
				{
				++m_reject_mega_fwd;
				continue;
				}

			if (target_mega_self_rev_score == FLT_MAX)
				{
				/////////////////////////////////////////////
				// db mega self-rev score is needed
				/////////////////////////////////////////////
				const uint8_t *target_mega_prof = target_data->m_mega_prof;
				const float *target_mega_pssm_rev = target_data->m_mega_pssm_rev;
				target_mega_self_rev_score = sw_flat_pssm_scoreonly(
					scratch_rows, scratch_pssms,
					target_mega_prof, LT, target_mega_pssm_rev, LT,
					m_params->m_feature_block_offsets,
					m_params->m_nfeat,
					-m_params->m_open,
					-m_params->m_ext);
				}

			/////////////////////////////////////////////
			// (reverse query) + db mega score
			/////////////////////////////////////////////
			float mega_rev_score = sw_flat_pssm_scoreonly(
				scratch_rows, scratch_pssms,
				target_mega_prof, LT, query_mega_pssm_rev, LQ,
				m_params->m_feature_block_offsets,
				m_params->m_nfeat,
				-m_params->m_open,
				-m_params->m_ext);

			/////////////////////////////////////////////
			// query self-rev mega from cache
			/////////////////////////////////////////////
			float query_mega_self_rev_score = m_query_mega_self_rev_scores[qidx];

			uint nmatch = path2posvecs3(path_buffer, ncol,
				lo_i, L_i, lo_j, L_j, pos_is, pos_js, m_maxL);

			/////////////////////////////////////////////
			// LDDT
			/////////////////////////////////////////////
			float lddt = flat_getlddt_muscle_some_floats(
				pos_is, L_i,
				pos_js, L_j,
				nmatch, distmx_i, distmx_j,
				considered_vec, preserved_vec);

			/////////////////////////////////////////////
			// DALI
			/////////////////////////////////////////////
			float dali = flat_get_dali4(
				pos_is, L_i,
				pos_js, L_j,
				nmatch, distmx_i, distmx_j);

			asserta(m_params->m_lddtx_w == 0);
			asserta(m_params->m_dalix_w == 0);

			asserta(target_mega_self_rev_score != FLT_MAX);
			float mega_self_score =
				(target_mega_self_rev_score + query_mega_self_rev_score)/2;
			
			/////////////////////////////////////////////
			// Test statistic (TS)
			/////////////////////////////////////////////
			float TS = 0;
			TS += mega_fwd_score;
			TS -= m_params->m_rev_w*mega_rev_score;
			TS -= m_params->m_self_w*mega_self_score;
			TS += m_params->m_nurev_w*nu_rev_score;	// TODO +ve sign?!
			TS += m_params->m_lddt_w*lddt*500;
			TS += m_params->m_dali_w*dali*10;

			if (TS < s_min_ts_fold)
				{
				++m_reject_min_ts;
				continue;
				}
			++m_accept_min_ts;

			hit.reset();
			hit.query = m_ptr_query_chains[qidx];
			hit.target = target_data->m_chain;
			hit.path = path_buffer;
			hit.ncol = ncol;
			hit.nu_fwd_score = float(nu_fwd_score);
			hit.nu_rev_score = float(nu_rev_score);
			hit.mega_fwd_score = mega_fwd_score;
			hit.mega_rev_score = mega_rev_score;
			hit.lddt = lddt;
			hit.dali = dali;
			hit.TS = TS;

			asserta(m_fhit);
			if (m_fhit)
				{
				// fprintf is thread-safe
				//fprintf(m_fhit, "%.3g\t%s\t%s\n",
				//	TS,
				//	query_label.c_str(),
				//	target_label.c_str());
				string str;
				str = query_label;
				str += "\t" + target_label;
				Psa(str, "\t%.3g", TS);
				Psa(str, "\t%.3g", float(nu_fwd_score));
				Psa(str, "\t%.3g", float(nu_rev_score));
				Psa(str, "\t%.3g", float(mega_fwd_score));
				Psa(str, "\t%.3g", float(mega_rev_score));
				Psa(str, "\t%.3g", lddt);
				Psa(str, "\t%.3g", dali);
				str += "\n";
				fputs(str.c_str(), m_fhit);
				}
			}

		struct_data::free_struct_data(target_data);
		}
	}
