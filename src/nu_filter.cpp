#include "myutils.h"
#include "nu_filter.h"
#include "flat_nu_aligner.h"
#include "flat_params.h"
#include "flat_helpers.h"

#define WRITE_QUERY_NU_SELF_REV_SCORES	0
#define WRITE_DB_NU_SELF_REV_SCORES		0
#define WRITE_TS_TERMS					1

const BCAData *nu_filter::m_dbbca;
const flat_params *nu_filter::m_params;
uint nu_filter::m_query_nchain = 0;
const vector<string> *nu_filter::m_ptr_query_labels = 0;
parasail_profile_t **nu_filter::m_query_parasail_profs = 0;
parasail_profile_t **nu_filter::m_query_parasail_prof_revs = 0;
const uint *nu_filter::m_query_lengths = 0;
int *nu_filter::m_query_self_rev_scores = 0;
const unordered_map<uint, vector<uint> > *nu_filter::m_dbidx_to_qidxs = 0;
const vector<uint> *nu_filter::m_dbidxs = 0;
atomic<uint> nu_filter::m_next;
const uint nu_filter::m_maxL = 4000;
atomic<uint> nu_filter::m_npair;
atomic<uint> nu_filter::m_reject_fwd;
atomic<uint> nu_filter::m_reject_cmb;
atomic<uint> nu_filter::m_npass;
atomic<uint> nu_filter::m_reject_mega_fwd;
atomic<uint> nu_filter::m_nhit;
const sid_t **nu_filter::m_query_distmxs;
const float **nu_filter::m_query_mega_pssms;
const float **nu_filter::m_query_mega_pssm_revs;
float *nu_filter::m_query_mega_self_rev_scores;
uint nu_filter::m_ndbidxs;
NF_MODE nu_filter::m_mode = NF_invalid;
vector<uint> nu_filter::m_qidxs_all;
FILE *nu_filter::m_fhits;
mutex nu_filter::m_hits_lock;

#if WRITE_TS_TERMS
static FILE *g_ftsv;
static mutex g_tsv_lock;
#endif

void nu_filter::set_query_data(
	const vector<string> &labels,
	parasail_profile_t **parasail_profs,
	parasail_profile_t **parasail_prof_revs,
	const float **query_mega_pssms,
	const float **query_mega_pssm_revs,
	const sid_t **query_distmxs,
	const uint *lengths,
	uint nchain)
	{
	m_ptr_query_labels = &labels;
	m_query_parasail_profs = parasail_profs;
	m_query_parasail_prof_revs = parasail_prof_revs;
	m_query_mega_pssms = query_mega_pssms;
	m_query_mega_pssm_revs = query_mega_pssm_revs;
	m_query_distmxs = query_distmxs;
	m_query_lengths = lengths;
	m_query_nchain = nchain;
	}

void nu_filter::set_query_mega_self_rev_scores(
	uint8_t **query_mega_profs)
	{
	asserta(m_query_mega_self_rev_scores == 0);
	asserta(m_query_nchain > 0);
	m_query_mega_self_rev_scores = myalloc(float, m_query_nchain);
	float *scratch_rows = myalloc(float, 2*m_maxL + 2);
	const float **scratch_pssms = myalloc(const float *, m_params->m_nfeat);

	Progress("Query Mega self-scores...");
	for (uint qidx = 0; qidx < m_query_nchain; ++qidx)
		{
		const uint LQ = m_query_lengths[qidx];
		const float *query_mega_pssm_rev = m_query_mega_pssm_revs[qidx];
		const uint8_t *query_mega_prof = query_mega_profs[qidx];
		float score = sw_flat_pssm_scoreonly(
			scratch_rows, scratch_pssms,
			query_mega_prof, LQ, query_mega_pssm_rev, LQ,
			m_params->m_feature_block_offsets,
			m_params->m_nfeat,
			-m_params->m_open,
			-m_params->m_ext);
		m_query_mega_self_rev_scores[qidx] = score;
		}
	myfree(scratch_rows);
	myfree(scratch_pssms);
	Progress(" done.\n");
	}

void nu_filter::set_query_self_rev_scores(
	uint8_t **query_codeseq_nus)
	{
	asserta(m_query_self_rev_scores == 0);
	asserta(m_query_nchain > 0);
	m_query_self_rev_scores = myalloc(int, m_query_nchain);
	const int open = flat_nu_aligner::m_open;
	const int ext = flat_nu_aligner::m_ext;
	uint workspace_bytes =
		parasail_nomalloc_sw_striped_profile_avx2_256_16_workspace_bytes(m_maxL);
	uint8_t *workspace = myalloc(uint8_t, workspace_bytes);

	Progress("Query Nu self-scores...");
#if WRITE_QUERY_NU_SELF_REV_SCORES
	FILE *ftmp = CreateStdioFile("query_self_rev_scores.tmp");
#endif
	for (uint qidx = 0; qidx < m_query_nchain; ++qidx)
		{
		const uint LQ = m_query_lengths[qidx];
		const uint8_t *query_codeseq_nu = query_codeseq_nus[qidx];
		parasail_profile_t *query_para_prof_rev = m_query_parasail_prof_revs[qidx];
		int nu_rev_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
			query_para_prof_rev, (const char *) query_codeseq_nu, LQ,
			open, ext, workspace, workspace_bytes);
		m_query_self_rev_scores[qidx] = nu_rev_score;
#if WRITE_QUERY_NU_SELF_REV_SCORES
		fprintf(ftmp, "%s\t%d\n", (*m_ptr_query_labels)[qidx].c_str(), nu_rev_score);
#endif
		}
#if WRITE_QUERY_NU_SELF_REV_SCORES
	CloseStdioFile(ftmp);
#endif
	Progress(" done.\n");
	myfree(workspace);
	}

void nu_filter::static_thread_body(uint threadidx)
	{
	const flat_params &params = *m_params;
	const BCAData &dbbca = *m_dbbca;
	//const vector<uint> &dbidxs = *m_dbidxs;
	//const unordered_map<uint, vector<uint> > dbidx_to_qidxs = *m_dbidx_to_qidxs;
	const int open = flat_nu_aligner::m_open;
	const int ext = flat_nu_aligner::m_ext;
	//const uint ndb = uint(dbidxs.size());
	uint workspace_bytes =
		parasail_nomalloc_sw_striped_profile_avx2_256_16_workspace_bytes(m_maxL);
	uint8_t *workspace = myalloc(uint8_t, workspace_bytes);
	chaq_vecs2 cv;
	chaq::alloc_chaq_vecs2(cv, m_maxL);
	//const uint ndbidxs = uint(m_dbidxs->size());
	const float revw = m_params->m_nu_filter_rev_w;
	const float selfw = m_params->m_nu_filter_self_w;
	float *scratch_rows = myalloc(float, 2*m_maxL + 2);
	const float **scratch_pssms = myalloc(const float *, m_params->m_nfeat);
	uint8_t *TB = myalloc(uint8_t, m_maxL*m_maxL);
	char *path_buffer = myalloc(char, 2*m_maxL);
	uint *pos_is = myalloc(uint, m_maxL);
	uint *pos_js = myalloc(uint, m_maxL);
	uint *considered_vec = myalloc(uint, m_maxL);
	uint *preserved_vec = myalloc(uint, m_maxL);
	uint scratch_buffer_bytes = 2*m_maxL;
	uint8_t *scratch_buffer = myalloc(uint8_t, scratch_buffer_bytes);
	chaq_vecs2 cv2;
	chaq::alloc_chaq_vecs2(cv2, m_maxL);

	for (;;)
		{
		uint k = m_next++;
		if (k + 1 == m_ndbidxs || (k%100 == 0 && k > 0 && k + 1 < m_ndbidxs))
			{
			static mutex progress_lock;
			progress_lock.lock();
			ProgressStep(k, m_ndbidxs, "Nu filter");
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
		struct_data *db_data = dbbca.get_struct_data(
			params, dbidx,
			&cv, scratch_buffer, scratch_buffer_bytes);
		const string &db_label = dbbca.m_Labels[dbidx];
		uint8_t *db_codeseq_nu = db_data->m_codeseq_nu;
		uint8_t *db_codeseq_nu_rev = db_data->m_codeseq_nu_rev;
		uint8_t *db_mega_prof = db_data->m_mega_prof;
		const sid_t *db_distmx = db_data->m_distmx;
		const uint LT = db_data->m_chain->get_length();

		parasail_profile_t *db_para_prof = db_data->m_parasail_prof;
		int db_self_rev_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
			db_para_prof, (const char *) db_codeseq_nu_rev, LT, open, ext,
			workspace, workspace_bytes);
#if WRITE_DB_NU_SELF_REV_SCORES
		static FILE *ftmp = 0;
		static mutex tmp_lock;
		tmp_lock.lock();
		if (ftmp == 0) ftmp = CreateStdioFile("db_nu_self_rev_scores.tmp");
		fprintf(ftmp, "%s\t%d\n", m_dbbca->m_Labels[dbidx].c_str(), db_self_rev_score);
		tmp_lock.unlock();
#endif

		float db_mega_self_rev_score = FLT_MAX; // calculate only if needed
		for (uint j = 0; j < nq; ++j)
			{
			++m_npair;
			uint qidx = qidxs[j];
			asserta(qidx < m_query_nchain);
			const string &query_label = (*m_ptr_query_labels)[qidx];
			parasail_profile_t *query_para_prof = m_query_parasail_profs[qidx];
			int fwd_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
				query_para_prof, (const char *) db_codeseq_nu, LT, open, ext,
				workspace, workspace_bytes);
			if (fwd_score < params.m_nu_filter_min_fwd_score)
				{
				++m_reject_fwd;
				continue;
				}

			parasail_profile_t *query_para_prof_rev = m_query_parasail_prof_revs[qidx];
			int nu_rev_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
				query_para_prof_rev, (const char *) db_codeseq_nu, LT, open, ext,
				workspace, workspace_bytes);

			float self_score = (db_self_rev_score + 
				m_query_self_rev_scores[qidx])/2.0f;

			float combined_score = float(fwd_score) - selfw*self_score - revw*nu_rev_score;
			if (combined_score < m_params->m_nu_filter_min_combined_score)
				{
				++m_reject_cmb;
				continue;
				}
			++m_npass;

			const uint LQ = m_query_lengths[qidx];
			const float *query_mega_pssm = m_query_mega_pssms[qidx];
			const float *query_mega_pssm_rev = m_query_mega_pssm_revs[qidx];

			uint lo_i, lo_j, ncol;
			float mega_fwd_score = sw_flat_pssm(
				scratch_rows, TB, scratch_pssms,
				db_mega_prof, LT, query_mega_pssm, LQ,
				m_params->m_feature_block_offsets,
				m_params->m_nfeat,
				-m_params->m_open,
				-m_params->m_ext,
				lo_i, lo_j, path_buffer, ncol);
			const uint L_i = LT;
			const uint L_j = LQ;
			const sid_t *distmx_i = db_distmx;
			const sid_t *distmx_j = m_query_distmxs[qidx];

			if (mega_fwd_score < m_params->m_mega_filter_min_fwd)
				{
				++m_reject_mega_fwd;
				continue;
				}

			if (db_mega_self_rev_score == FLT_MAX)
				{
				const uint8_t *db_mega_prof = db_data->m_mega_prof;
				const float *db_mega_pssm_rev = db_data->m_mega_pssm_rev;
				db_mega_self_rev_score = sw_flat_pssm_scoreonly(
					scratch_rows, scratch_pssms,
					db_mega_prof, LT, db_mega_pssm_rev, LT,
					m_params->m_feature_block_offsets,
					m_params->m_nfeat,
					-m_params->m_open,
					-m_params->m_ext);
				}

			float mega_rev_score = sw_flat_pssm_scoreonly(
				scratch_rows, scratch_pssms,
				db_mega_prof, LT, query_mega_pssm_rev, LQ,
				m_params->m_feature_block_offsets,
				m_params->m_nfeat,
				-m_params->m_open,
				-m_params->m_ext);

			uint nmatch = path2posvecs3(path_buffer, ncol,
				lo_i, L_i, lo_j, L_j, pos_is, pos_js, m_maxL);

			float lddt = flat_getlddt_muscle_some_floats(
				pos_is, L_i,
				pos_js, L_j,
				nmatch, distmx_i, distmx_j,
				considered_vec, preserved_vec);

			float dali = flat_get_dali4(
				pos_is, L_i,
				pos_js, L_j,
				nmatch, distmx_i, distmx_j);

			asserta(m_params->m_lddtx_w == 0);
			asserta(m_params->m_dalix_w == 0);

			float query_mega_self_rev_score = m_query_mega_self_rev_scores[qidx];
			asserta(db_mega_self_rev_score != FLT_MAX);
			float mega_self_score =
				(db_mega_self_rev_score + query_mega_self_rev_score)/2;
			
			float TS = 0;
			TS += mega_fwd_score;
			TS -= m_params->m_rev_w*mega_rev_score;
			TS -= m_params->m_self_w*mega_self_score;
			TS += m_params->m_nurev_w*nu_rev_score;	// TODO +ve sign?!
			TS += m_params->m_lddt_w*lddt*500;
			TS += m_params->m_dali_w*dali*10;

			if (m_fhits != 0)
				{
				m_hits_lock.lock();
				fprintf(m_fhits, "%.3g\t%s\t%s\n",
					TS,
					query_label.c_str(),
					db_label.c_str());
				m_hits_lock.unlock();
				}

#if WRITE_TS_TERMS
			g_tsv_lock.lock();
			fprintf(g_ftsv, "%s\t%s", query_label.c_str(), db_label.c_str());
			fprintf(g_ftsv, "\t%.3g", mega_fwd_score);
			fprintf(g_ftsv, "\t%.3g", mega_rev_score);
			fprintf(g_ftsv, "\t%d", nu_rev_score);
			fprintf(g_ftsv, "\t%.3g", mega_self_score);
			fprintf(g_ftsv, "\t%.3g", lddt);
			fprintf(g_ftsv, "\t%.3g", dali);
			fprintf(g_ftsv, "\t%.3g", TS);
			fprintf(g_ftsv, "\n");
			g_tsv_lock.unlock();
#endif
			}

		struct_data::free_struct_data(db_data);
		}
	}

void nu_filter::run_filter()
	{
	asserta(m_params != 0);

	m_npair = 0;
	m_reject_fwd = 0;
	m_reject_cmb = 0;
	m_npass = 0;
	m_next = 0;

	ProgressStep(0, m_ndbidxs, "Nu filter");
#if WRITE_TS_TERMS
	g_ftsv = CreateStdioFile("nu_filter.tmp");
#endif

	vector<thread *> ts;
	uint ThreadCount = GetRequestedThreadCount();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(static_thread_body, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	ProgressLog("%10u  Nu filter npair\n", m_npair.load());
	ProgressLog("%10u  Nu filter nreject_fwd\n", m_reject_fwd.load());
	ProgressLog("%10u  Nu filter nreject_cmb\n", m_reject_cmb.load());
	ProgressLog("%10u  Nu filter pass\n", m_npass.load());
#if WRITE_TS_TERMS
	CloseStdioFile(g_ftsv);
#endif
	}

void nu_filter::run_filter_all_vs_all(const BCAData &dbbca)
	{
	m_ndbidxs = dbbca.GetChainCount();
	m_dbbca = &dbbca;
	m_dbidxs = 0;
	m_dbidx_to_qidxs = 0;
	m_mode = NF_all_vs_all;
	m_qidxs_all.clear();
	m_qidxs_all.reserve(m_query_nchain);
	for (uint i = 0; i < m_query_nchain; ++i)
		m_qidxs_all.push_back(i);

	run_filter();
	}

void nu_filter::run_filter_post_kappa(
	const BCAData &dbbca,
	const vector<uint> &dbidxs,
	const unordered_map<uint, vector<uint> > &dbidx_to_qidxs)
	{
	m_ndbidxs = uint(dbidxs.size());
	m_dbbca = &dbbca;
	m_dbidxs = &dbidxs;
	m_dbidx_to_qidxs = &dbidx_to_qidxs;
	m_mode = NF_kappa;
	m_qidxs_all.clear();

	run_filter();
	}
