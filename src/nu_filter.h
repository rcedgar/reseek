#pragma once

#include "parasail_nomalloc.h"
#include "bcadata.h"

class nu_filter
	{
public:
	static const uint m_maxL;
	static const flat_params *m_params;
	static uint m_query_nchain;
	static const vector<string> *m_ptr_query_labels;
	static parasail_profile_t **m_query_parasail_profs;
	static parasail_profile_t **m_query_parasail_prof_revs;
	static const uint *m_query_lengths;
	static int *m_query_self_rev_scores;
	static atomic<uint> m_next;
	static const BCAData *m_dbbca;
	static const unordered_map<uint, vector<uint> > *m_dbidx_to_qidxs;
	static const vector<uint> *m_dbidxs;

	static atomic<uint> m_npair;
	static atomic<uint> m_reject_fwd;
	static atomic<uint> m_reject_cmb;
	static atomic<uint> m_npass;

public:
	static void set_query_parasail_profiles(
		const vector<string> &labels,
		parasail_profile_t **query_parasail_profs,
		parasail_profile_t **query_parasail_prof_revs,
		const uint *lengths,
		uint nchain);

	static void set_query_self_rev_scores(
		uint8_t **query_codeseq_nus);

	static void run_filter(
		const flat_params &params,
		const BCAData &dbbca,
		const vector<uint> &dbidxs,
		const unordered_map<uint, vector<uint> > &dbidx_to_qidxs);

	static void static_thread_body(uint threadidx);
	};
