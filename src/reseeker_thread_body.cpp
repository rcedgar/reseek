#include "myutils.h"
#include "reseeker.h"
#include "bcadata_struct.h"
#include "flat_nu_aligner.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "hitdata.h"
#include "sw_flat_pssm_xdrop.h"
#include "rankedscoresbag.h"
#include "kappa_hsp.h"

static void validate_hit(const hitdata &hit)
	{
	const char *qaa = hit.query->m_aa->m_data;
	const char *taa = hit.target->m_aa->m_data;
	uint qpos = hit.qlo;
	uint tpos = hit.tlo;
	const uint LQ = hit.query->m_L;
	const uint LT = hit.target->m_L;
	for (uint i = 0; i < hit.ncol; ++i)
		{
		char c = hit.path[i];
		switch (c)
			{
		case 'M':
			asserta(qpos < LQ);
			asserta(tpos < LT);
			++qpos;
			++tpos;
			break;

		case 'D':
			asserta(tpos < LT);
			++tpos;
			break;

		case 'I':
			asserta(qpos < LQ);
			++qpos;
			break;

		default:
			asserta(false);
			}
		}

	}

void reseeker::static_thread_body(uint threadidx)
	{
	const BCAData &dbbca = *m_dbbca;

	const int nu_open = flat_nu_aligner::m_open;
	const int nu_ext = flat_nu_aligner::m_ext;

	const flat_params &params = *m_params;
	const float revw = m_params->m_nu_filter_rev_w;
	const float selfw = m_params->m_nu_filter_self_w;
	const bool nu_only = m_params->m_nu_only;

	uint workspace_bytes =
		parasail_nomalloc_sw_striped_profile_avx2_256_16_workspace_bytes(flat_params::m_maxL);
	uint scratch_buffer_bytes = 2*flat_params::m_maxL;

	uint8_t *workspace = myalloca(uint8_t, workspace_bytes);
	float *scratch_rows = myalloc(float, 2*flat_params::m_maxL + 3);
	const float **scratch_pssms = myalloc(const float *, m_params->m_nfeat);
	uint8_t *TB = myalloc(uint8_t, flat_params::m_maxL*flat_params::m_maxL);
	uint8_t *TB_bwd = myalloc(uint8_t, flat_params::m_maxL*flat_params::m_maxL);
	// 3 regions for sw_flat_pssm_xdrop_hsp (fwd/bwd/merged)
	char *path_buffer = myalloc(char, 6*flat_params::m_maxL + 16);
	char *path_check = myalloc(char, 2*flat_params::m_maxL);
	uint *pos_is = myalloc(uint, flat_params::m_maxL);
	uint *pos_js = myalloc(uint, flat_params::m_maxL);
	uint *considered_vec = myalloc(uint, flat_params::m_maxL);
	uint *preserved_vec = myalloc(uint, flat_params::m_maxL);
	uint8_t *scratch_buffer = myalloc(uint8_t, scratch_buffer_bytes);
	chaq_vecs2 cv;
	chaq::alloc_chaq_vecs2(cv, flat_params::m_maxL);
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
		const vector<uint> *ptr_diagscores = 0;
		const vector<uint> *ptr_hspdiags = 0;
		const vector<uint> *ptr_hsplos = 0;
		const vector<uint> *ptr_hsplens = 0;
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
			unordered_map<uint, vector<uint> >::const_iterator iter_diag =
				m_dbidx_to_diagscores->find(dbidx);
			asserta(iter_diag != m_dbidx_to_diagscores->end());
			ptr_diagscores = &iter_diag->second;
			asserta(ptr_diagscores->size() == ptr_qidxs->size());
			if (m_dbidx_to_hspdiags != 0)
				{
				asserta(m_dbidx_to_hsplos != 0 && m_dbidx_to_hsplens != 0);
				auto it_d = m_dbidx_to_hspdiags->find(dbidx);
				auto it_lo = m_dbidx_to_hsplos->find(dbidx);
				auto it_ln = m_dbidx_to_hsplens->find(dbidx);
				asserta(it_d != m_dbidx_to_hspdiags->end());
				asserta(it_lo != m_dbidx_to_hsplos->end());
				asserta(it_ln != m_dbidx_to_hsplens->end());
				ptr_hspdiags = &it_d->second;
				ptr_hsplos = &it_lo->second;
				ptr_hsplens = &it_ln->second;
				asserta(ptr_hspdiags->size() == ptr_qidxs->size());
				}
			break;
			}

		default: asserta(false);
			}

		const vector<uint> &qidxs = *ptr_qidxs;
		const uint nq = uint(qidxs.size());
		asserta(nq > 0);
		struct_data *target_data = bca_get_struct_data(dbbca,
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
			uint kappa_diag_score = 0;
			if (m_mode == NF_kappa)
				kappa_diag_score = (*ptr_diagscores)[j];
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
			float nu_combined_score = float(nu_fwd_score) - selfw*self_score - revw*nu_rev_score;
			if (nu_combined_score < m_params->m_nu_filter_min_combined_score)
				{
				++m_nu_reject_cmb;
				continue;
				}
			++m_npass;

			if (nu_only)
				{
				if (m_fhit)
					{
					string str;
					str = query_label;
					str += "\t" + target_label;
					Psa(str, "\t%.3g", nu_combined_score);
					str += "\n";
					// fprintf & fputs are thread-safe
					fputs(str.c_str(), m_fhit);
					}
				continue;
				}

			const uint LQ = m_query_lengths[qidx];
			const float *query_mega_pssm = m_query_mega_pssms[qidx];
			const float *query_mega_pssm_rev = m_query_mega_pssm_revs[qidx];

			/////////////////////////////////////////////
			// Mega forward score
			/////////////////////////////////////////////
			uint lo_i, lo_j, ncol;
			float mega_fwd_score = 0;
			bool used_hsp_align = false;

			const uint gate = flat_params::m_hsp_align_min_length;
			const bool long_pair = (LQ >= gate || LT >= gate);
			uint hsp_diag = HSP_SEED_NONE;
			uint hsp_lo = 0;
			uint hsp_len = 0;
			if (ptr_hspdiags != 0)
				{
				hsp_diag = (*ptr_hspdiags)[j];
				hsp_lo = (*ptr_hsplos)[j];
				hsp_len = (*ptr_hsplens)[j];
				}
			const bool have_seed = (hsp_diag != HSP_SEED_NONE && hsp_len > 0);
			const bool try_hsp = flat_params::m_hsp_align && long_pair && have_seed;
			if (flat_params::m_hsp_align && long_pair && !have_seed)
				++m_hsp_align_no_seed;

			if (try_hsp)
				{
				++m_hsp_align_try;
				int mini = 0, minj = 0, n = 0;
				kappa_get_hsp_limits(int(LQ), int(LT), int(hsp_diag),
					mini, minj, n);
				(void) n;
				// Seed mid-HSP (Mu XDropHSP style), not the HSP start edge.
				const int seed_off = int(hsp_lo) + int(hsp_len) / 2;
				const int pos_query = mini + seed_off;
				const int pos_target = minj + seed_off;
				asserta(pos_query >= 0 && uint(pos_query) < LQ);
				asserta(pos_target >= 0 && uint(pos_target) < LT);

				// sw_flat_pssm: prof = target, pssm = query
				mega_fwd_score = sw_flat_pssm_xdrop_hsp(
					scratch_rows, TB, TB_bwd, scratch_pssms,
					target_mega_prof, LT, query_mega_pssm, LQ,
					m_params->m_feature_block_offsets,
					m_params->m_nfeat,
					uint(pos_target), uint(pos_query),
					flat_params::m_hsp_x2,
					-m_params->m_open,
					-m_params->m_ext,
					lo_i, lo_j, path_buffer, ncol);

				if (flat_params::m_hsp_align_check)
					{
					uint lo_i_f, lo_j_f, ncol_f;
					float full_score = sw_flat_pssm(
						scratch_rows, TB_bwd, scratch_pssms,
						target_mega_prof, LT, query_mega_pssm, LQ,
						m_params->m_feature_block_offsets,
						m_params->m_nfeat,
						-m_params->m_open,
						-m_params->m_ext,
						lo_i_f, lo_j_f, path_check, ncol_f);
					(void) ncol_f;
					const float eps = 1e-3f + 1e-5f * max(fabsf(full_score), 1.0f);
					if (mega_fwd_score > full_score + eps)
						{
						++m_hsp_align_check_bug;
						Die("hsp_align_check: xdrop=%.6g > full=%.6g (q=%s t=%s diag=%u lo=%u len=%u)",
							mega_fwd_score, full_score,
							query_label.c_str(), target_label.c_str(),
							hsp_diag, hsp_lo, hsp_len);
						}
					else if (fabsf(mega_fwd_score - full_score) <= eps)
						++m_hsp_align_check_ok;
					else
						++m_hsp_align_check_xdrop_lt;
					}

				if (mega_fwd_score < m_params->m_mega_filter_min_fwd || ncol == 0)
					{
					// Cheap reject — do not fall back to full Mega SW.
					++m_hsp_align_fallback;
					++m_hsp_align_reject_score;
					++m_reject_mega_fwd;
					continue;
					}
				used_hsp_align = true;
				}
			else
				{
				mega_fwd_score = sw_flat_pssm(
					scratch_rows, TB, scratch_pssms,
					target_mega_prof, LT, query_mega_pssm, LQ,
					m_params->m_feature_block_offsets,
					m_params->m_nfeat,
					-m_params->m_open,
					-m_params->m_ext,
					lo_i, lo_j, path_buffer, ncol);
				}
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
				lo_i, L_i, lo_j, L_j, pos_is, pos_js, flat_params::m_maxL);
			if (nmatch == 0)
				{
				++m_reject_mega_fwd;
				if (used_hsp_align)
					{
					++m_hsp_align_fallback;
					++m_hsp_align_reject_path;
					}
				continue;
				}
			if (used_hsp_align)
				++m_hsp_align_used;

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

			if (TS < m_mints)
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
			hit.qlo = lo_j;
			hit.tlo = lo_i;
			hit.TS = TS;
			hit.fill(params);
#if DEBUG
			validate_hit(hit);
#endif
			write_hit(params, hit);
			}

		struct_data::free_struct_data(target_data);
		}
	}
