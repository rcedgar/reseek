#pragma once

#include "parasail_nomalloc.h"
#include "bcadata.h"

enum NF_MODE
	{
	NF_invalid,
	NF_all_vs_all,
	NF_kappa
	};

class reseeker
	{
public:
	static const uint m_maxL;
	static const flat_params *m_params;
	static const flat_params *m_params2;
	static const flat_params *m_params3;
	static uint m_query_nchain;
	static const vector<string> *m_ptr_query_labels;
	static parasail_profile_t **m_query_parasail_profs;
	static parasail_profile_t **m_query_parasail_prof_revs;
	static const sid_t **m_query_distmxs;
	static const uint *m_query_lengths;
	static const float **m_query_mega_pssms;
	static const float **m_query_mega_pssm_revs;
	static int *m_query_self_rev_scores;
	static float *m_query_mega_self_rev_scores;
	static atomic<uint> m_next;
	static const BCAData *m_dbbca;
	static const unordered_map<uint, vector<uint> > *m_dbidx_to_qidxs;
	static const vector<uint> *m_dbidxs;
	static uint m_ndbidxs;
	static NF_MODE m_mode;
	static vector<uint> m_qidxs_all;

	static atomic<uint> m_npair;
	static atomic<uint> m_reject_fwd;
	static atomic<uint> m_reject_cmb;
	static atomic<uint> m_npass;
	static atomic<uint> m_reject_mega_fwd;
	static atomic<uint> m_reject_min_fold_ts;
	static atomic<uint> m_accept_min_fold_ts;
	static atomic<uint> m_nhit;

	static FILE *m_fhits;
	static mutex m_hits_lock;

public:
	static void set_params(const flat_params &params)
		{
		m_params = &params;
		}

	static void set_query_data(
		const vector<string> &labels,
		parasail_profile_t **query_parasail_profs,
		parasail_profile_t **query_parasail_prof_revs,
		const float **query_mega_pssms,
		const float **query_mega_pssm_revs,
		const sid_t **query_distmxs,
		const uint *query_lengths,
		uint nchain);

	static void set_query_self_rev_scores(
		uint8_t **query_codeseq_nus);

	static void set_query_mega_self_rev_scores(
		uint8_t **query_mega_profs);

	static void run_filter();

	static void run_filter_all_vs_all(
		const BCAData &dbbca);

	static void run_filter_post_kappa(
		const BCAData &dbbca,
		const vector<uint> &dbidxs,
		const unordered_map<uint, vector<uint> > &dbidx_to_qidxs);

	static void static_thread_body(uint threadidx);
	};
