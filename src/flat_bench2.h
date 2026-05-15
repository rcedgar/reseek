#pragma once

#include "fastbench.h"
#include "chain_data.h"
#include "flat_alphas.h"

struct flat_bench2_thread_data
	{
	const float **m_scratch_pssms = 0;
	float *m_scratch_rows = 0;
	uint8_t *m_TB = 0;
	char *m_path_buffer = 0;
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
	atomic<uint> m_aligned_pair_count = 0;

public:
	void search(uint nthread, bool pin_threads);
	void align_pair(uint pairidx, flat_bench2_thread_data &TD);
	void load_chains(const vector<flat_chain_t *> &chains);
	void thread_body(uint ThreadIdx);
	void update_params(
		const vector<string> &names,
		const vector<float> &values);

public:
	static void static_thread_body(flat_bench2 *SB, uint threadidx);
	};
