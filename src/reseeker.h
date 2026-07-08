#pragma once

#include "parasail_nomalloc.h"
#include "bcadata.h"
#include "userfields.h"

class flat_chain_t;
class hitdata;

enum NF_MODE
	{
	NF_invalid,
	NF_all_vs_all,
	NF_kappa
	};

typedef void (*ptr_thread_body_fn)(uint threadidx);

class reseeker
	{
public:
	static const flat_params *m_params;
	static uint m_query_nchain;
	static const flat_chain_t **m_ptr_query_chains;
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
	static const unordered_map<uint, vector<uint> > *m_dbidx_to_diagscores;
	static const vector<uint> *m_dbidxs;
	static uint m_ndbidxs;
	static NF_MODE m_mode;
	static vector<uint> m_qidxs_all;
	static float m_mints;
	static uint m_max_queries_per_target;

	static atomic<uint> m_npair;
	static atomic<uint> m_reject_fwd;
	static atomic<uint> m_nu_reject_cmb;
	static atomic<uint> m_npass;
	static atomic<uint> m_reject_mega_fwd;
	static atomic<uint> m_reject_min_ts;
	static atomic<uint> m_accept_min_ts;
	static atomic<uint> m_nhit;

	static vector<USERFIELD> m_UFs;

	static FILE *m_fhit;
	static FILE *m_faln;
	//static mutex m_hit_lock; // exploit fputs thread-safety
	static mutex m_aln_lock;

public:
	static void set_params(const flat_params &params)
		{
		m_params = &params;
		init_userfields();
		}

	static void set_query_data(
		const flat_chain_t **ptr_query_chains,
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

	static void search();

	static void search_all_vs_all(
		const BCAData &dbbca);

	static void search_post_kappa(
		const BCAData &dbbca,
		const vector<uint> &dbidxs,
		const unordered_map<uint, vector<uint> > &dbidx_to_qidxs,
		const unordered_map<uint, vector<uint> > &dbidx_to_diagscores);

	static void static_thread_body(uint threadidx);
	//static void static_thread_body_nusort(uint threadidx);
	static void write_hit(const flat_params &params, const hitdata &hit);
	static void init_userfields();
	static void append_userfield(
		string &s,
		const hitdata &hit,
		USERFIELD UF);

	static void close_files()
		{
		CloseStdioFile(m_fhit);
		CloseStdioFile(m_faln);
		}

private:
	static void write_tsv(const hitdata &hit);
	static void write_aln(const hitdata &hit);
	};
