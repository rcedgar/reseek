#pragma once

#include "fastbench.h"
#include "chain_data.h"
#include "flat_alphas.h"
#include "parasail.h"

class flat_bench2_thread_data
	{
public:
	const float **m_scratch_pssms = 0;
	float *m_scratch_rows = 0;
	uint8_t *m_TB = 0;
	char *m_path_buffer = 0;
	float *m_colscores = 0;
	parasail_result_t *m_parasail_result = 0;

private:
	flat_bench2_thread_data() = delete;

public:
	flat_bench2_thread_data(uint maxL, uint nfeat)
		{
		m_scratch_rows = myalloc(float, 2*maxL + 2);
		m_scratch_pssms = myalloc(const float *, nfeat);
		m_TB = myalloc(uint8_t, maxL*maxL);
		m_path_buffer = myalloc(char, 2*maxL);
		m_colscores = myalloc(float, 2*maxL);
		m_parasail_result = 0;
		}

	~flat_bench2_thread_data()
		{
		myfree(m_scratch_pssms);
		myfree(m_scratch_rows);
		myfree(m_TB);
		myfree(m_path_buffer);
		if (m_parasail_result != 0)
			parasail_result_free(m_parasail_result);
		}
	};

class flat_bench2 : public FastBench
	{
public:
// Hard-coded score-like not Evalue-like
//	SBSCORE m_SBS = SBS_Evalue;

	static chain_data **m_cdvec;
	static uint m_maxL;

public:
	atomic<uint> m_next_pairidx = 0;
	float *m_self_rev_scores = 0;
	float *m_nu_self_rev_scores = 0;
	bool m_nu_filter = true;
	bool m_nu_only = false;
	atomic<uint> m_aln_count = 0;
	atomic<uint> m_mega_fwd_reject_count= 0;
	atomic<uint> m_mu_fwd_reject_count = 0;
	atomic<uint> m_mu_combined_reject_count= 0;

public:
	void search(uint nthread, bool pin_threads);
	void align_pair(uint pairidx, flat_bench2_thread_data &TD);
	void load_chains(const vector<flat_chain_t *> &chains);
	void thread_body(uint ThreadIdx);
	void update_params(
		const vector<string> &names,
		const vector<float> &values);
	void set_self_rev_scores();
	void set_nu_self_rev_scores();

public:
	static void static_thread_body(flat_bench2 *SB, uint threadidx);
	};
