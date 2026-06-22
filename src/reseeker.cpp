#include "myutils.h"
#include "reseeker.h"
#include "flat_nu_aligner.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "hitdata.h"
#include "reseek_hit_sink.h"

uint reseeker::m_query_nchain = 0;
const BCAData *reseeker::m_dbbca = 0;
const flat_params *reseeker::m_params = 0;
const flat_chain_t **reseeker::m_ptr_query_chains = 0;
const vector<string> *reseeker::m_ptr_query_labels = 0;
parasail_profile_t **reseeker::m_query_parasail_profs = 0;
parasail_profile_t **reseeker::m_query_parasail_prof_revs = 0;
const uint *reseeker::m_query_lengths = 0;
chain_slice *reseeker::m_query_slices = 0;
int *reseeker::m_query_self_rev_scores = 0;
const unordered_map<uint, vector<uint> > *reseeker::m_dbidx_to_qidxs = 0;
const unordered_map<uint, vector<uint> > *reseeker::m_dbidx_to_diagscores = 0;
const vector<uint> *reseeker::m_dbidxs = 0;
atomic<uint> reseeker::m_next;
atomic<uint> reseeker::m_npair;
atomic<uint> reseeker::m_reject_fwd;
atomic<uint> reseeker::m_nu_reject_cmb;
atomic<uint> reseeker::m_npass;
atomic<uint> reseeker::m_reject_mega_fwd;
atomic<uint> reseeker::m_accept_min_ts;
atomic<uint> reseeker::m_reject_min_ts;
atomic<uint> reseeker::m_nhit;
const sid_t **reseeker::m_query_distmxs;
const float **reseeker::m_query_mega_pssms;
const float **reseeker::m_query_mega_pssm_revs;
float *reseeker::m_query_mega_self_rev_scores;
uint reseeker::m_ndbidxs;
NF_MODE reseeker::m_mode = NF_invalid;
vector<uint> reseeker::m_qidxs_all;
FILE *reseeker::m_fhit;
FILE *reseeker::m_faln;
//mutex reseeker::m_hit_lock; // exploit fputs thread-safety
mutex reseeker::m_aln_lock;
float reseeker::m_mints = 0;
uint reseeker::m_max_queries_per_target = 0;

void reseeker::set_query_data(
	const flat_chain_t **ptr_query_chains,
	const vector<string> &labels,
	parasail_profile_t **parasail_profs,
	parasail_profile_t **parasail_prof_revs,
	const float **query_mega_pssms,
	const float **query_mega_pssm_revs,
	const sid_t **query_distmxs,
	const uint *lengths,
	uint nchain)
	{
	m_ptr_query_chains = ptr_query_chains;
	m_ptr_query_labels = &labels;
	m_query_parasail_profs = parasail_profs;
	m_query_parasail_prof_revs = parasail_prof_revs;
	m_query_mega_pssms = query_mega_pssms;
	m_query_mega_pssm_revs = query_mega_pssm_revs;
	m_query_distmxs = query_distmxs;
	m_query_lengths = lengths;
	m_query_nchain = nchain;

	myfree(m_query_slices);
	m_query_slices = myalloc(chain_slice, nchain);
	for (uint i = 0; i < nchain; ++i)
		m_query_slices[i] = chain_slice_identity(i, lengths[i]);
	}

void reseeker::set_query_mega_self_rev_scores(
	uint8_t **query_mega_profs)
	{
	if (m_params->m_nu_only) return;
	asserta(m_query_mega_self_rev_scores == 0);
	asserta(m_query_nchain > 0);
	m_query_mega_self_rev_scores = myalloc(float, m_query_nchain);
	float *scratch_rows = myalloc(float, 2*flat_params::m_maxL + 2);
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

void reseeker::set_query_self_rev_scores(
	uint8_t **query_codeseq_nus)
	{
	asserta(m_query_self_rev_scores == 0);
	asserta(m_query_nchain > 0);
	m_query_self_rev_scores = myalloc(int, m_query_nchain);
	const int open = flat_nu_aligner::m_open;
	const int ext = flat_nu_aligner::m_ext;
	uint workspace_bytes =
		parasail_nomalloc_sw_striped_profile_avx2_256_16_workspace_bytes(flat_params::m_maxL);
	uint8_t *workspace = myalloc(uint8_t, workspace_bytes);

	Progress("Query Nu self-scores...");
	for (uint qidx = 0; qidx < m_query_nchain; ++qidx)
		{
		const uint LQ = m_query_lengths[qidx];
		const uint8_t *query_codeseq_nu = query_codeseq_nus[qidx];
		parasail_profile_t *query_para_prof_rev = m_query_parasail_prof_revs[qidx];
		int nu_rev_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
			query_para_prof_rev, (const char *) query_codeseq_nu, LQ,
			open, ext, workspace, workspace_bytes);
		m_query_self_rev_scores[qidx] = nu_rev_score;
		}
	Progress(" done.\n");
	myfree(workspace);
	}

void reseeker::search()
	{
	asserta(m_params != 0);

	m_npair = 0;
	m_reject_fwd = 0;
	m_nu_reject_cmb = 0;
	m_npass = 0;
	m_next = 0;

	vector<FILE *> fs;
	if (optset_output)
		reseeker::m_fhit = CreateStdioFile(opt(output));

	const uint ThreadCount = GetRequestedThreadCount();
	reseek_hit_sink_begin(ThreadCount);

	if (optset_max_nu_accepts)
		flat_params::m_max_nu_filter_accepts = opt(max_nu_accepts);
	else
		flat_params::m_max_nu_filter_accepts = 0;
	ptr_thread_body_fn thread_body =
		(flat_params::m_max_nu_filter_accepts > 0 ?
		static_thread_body_nusort :
		static_thread_body);

	ProgressStep(0, m_ndbidxs, "reseek");

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		//thread *t = new thread(static_thread_body, ThreadIndex);
		//thread *t = new thread(static_thread_body_nusort, ThreadIndex);
		thread *t = new thread(thread_body, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	reseek_hit_sink_flush(reseeker::m_fhit);
	reseeker::close_files();

	ProgressLog("%10u  Nu filter max accepts\n",
		flat_params::m_max_nu_filter_accepts);
	ProgressLog("%10u  Nu filter npair\n", m_npair.load());
	ProgressLog("%10u  Nu filter nreject_fwd\n", m_reject_fwd.load());
	ProgressLog("%10u  Nu filter nreject_cmb\n", m_nu_reject_cmb.load());
	ProgressLog("%10u  Nu filter pass\n", m_npass.load());
	ProgressLog("%10u  Mega filter reject\n", m_reject_min_ts.load());
	ProgressLog("%10u  Mega filter pass\n", m_accept_min_ts.load());
	}

void reseeker::search_all_vs_all(const BCAData &dbbca)
	{
	m_ndbidxs = dbbca.GetChainCount();
	m_dbbca = &dbbca;
	m_dbidxs = 0;
	m_dbidx_to_qidxs = 0;
	m_dbidx_to_diagscores = 0;
	m_mode = NF_all_vs_all;
	m_qidxs_all.clear();
	m_qidxs_all.reserve(m_query_nchain);
	for (uint i = 0; i < m_query_nchain; ++i)
		m_qidxs_all.push_back(i);

	search();
	}

void reseeker::search_post_kappa(
	const BCAData &dbbca,
	const vector<uint> &dbidxs,
	const unordered_map<uint, vector<uint> > &dbidx_to_qidxs,
	const unordered_map<uint, vector<uint> > &dbidx_to_diagscores)
	{
	m_ndbidxs = uint(dbidxs.size());
	m_dbbca = &dbbca;
	m_dbidxs = &dbidxs;
	m_dbidx_to_qidxs = &dbidx_to_qidxs;
	m_dbidx_to_diagscores = &dbidx_to_diagscores;
	m_mode = NF_kappa;
	m_qidxs_all.clear();

	search();
	}
