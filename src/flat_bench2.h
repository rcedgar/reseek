#pragma once

#include "fastbench.h"
#include "chain_data.h"
#include "flat_alphas.h"
#include "parasail.h"
#include "flat_bench2_thread_data.h"

class flat_bench2 : public FastBench
	{
public:
// Hard-coded score-like not Evalue-like
//	SBSCORE m_SBS = SBS_Evalue;

	static uint m_maxL;

public:
	chain_data **m_cdvec = 0;
	atomic<uint> m_next_pairidx = 0;
	atomic<uint> m_next_domidx = 0;
	float *m_self_rev_scores = 0;
	float *m_nu_self_rev_scores = 0;
	bool m_nu_only = false;
	bool m_timealn = false;
	atomic<uint> m_aln_count = 0;
	atomic<uint> m_mega_fwd_test_count= 0;
	atomic<uint> m_mega_fwd_pass_count= 0;
	atomic<uint> m_mu_fwd_reject_count = 0;
	atomic<uint> m_mu_combined_reject_count = 0;

	bool m_output_nu_paths = false;
	bool m_input_mega_paths = false;

	vector<uint> m_mega_path_is;
	vector<uint> m_mega_path_js;
	vector<uint> m_mega_path_lo_i_fwds;
	vector<uint> m_mega_path_lo_j_fwds;
	vector<uint> m_mega_path_lo_i_revs;
	vector<uint> m_mega_path_lo_j_revs;
	vector<string> m_mega_path_cigar_fwds;
	vector<string> m_mega_path_cigar_revs;
	vector<float> m_mega_path_score_fwds;
	vector<float> m_mega_path_score_revs;

public:
	void search(uint nthread, bool pin_threads);
	void align_pair(uint pairidx, flat_bench2_thread_data &TD);
	void align_pair_timealn(uint pairidx, flat_bench2_thread_data &TD);
	void align_pair_nu_only(uint pairidx, flat_bench2_thread_data &TD);
	void align_pair_output_nu_paths(uint pairidx, flat_bench2_thread_data &TD);
	void align_pair_input_mega_paths(uint pairidx, flat_bench2_thread_data &TD);
	void load_chains(const vector<flat_chain_t *> &chains);
	void thread_body(uint ThreadIdx);
	void thread_body_set_mega_self_rev_scores(uint ThreadIdx);
	void update_params(
		const vector<string> &names,
		const vector<float> &values);
	void set_mega_self_rev_score(
		uint domidx, flat_bench2_thread_data &TD);
	void set_mega_self_rev_scores();
	void set_nu_self_rev_scores();
	void load_mega_paths(const string &fn);

public:
	static void static_thread_body(
		flat_bench2 *SB, uint threadidx);

	static void static_thread_body_set_mega_self_rev_scores(
		flat_bench2 *SB, uint threadidx);

	float calc_ts(
		uint i, uint j,
		const chain_data &cd_i,
		const chain_data &cd_j,
		uint fwd_lo_i, uint fwd_lo_j,
		const char *fwd_path,
		uint fwd_ncol,
		uint rev_lo_i, uint rev_lo_j,
		const char *rev_path,
		uint rev_ncol,
		flat_bench2_thread_data &TD);

	static float score_path(
		const chain_data &cd_i,
		uint lo_i,
		const chain_data &cd_j,
		uint lo_j,
		const char *path,
		uint ncol);

	static float score_pos_pair(
		const uint8_t *mega_prof_i, uint pos_i, uint L_i,
		const uint8_t *mega_prof_j, uint pos_j, uint L_j);
	};
