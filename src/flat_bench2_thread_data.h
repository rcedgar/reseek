#pragma once
#include "flat_params.h"
#include "parasail_nomalloc.h"

class flat_bench2_thread_data
	{
public:
	const float **m_scratch_pssms = 0;
	float *m_scratch_rows = 0;
	uint8_t *m_TB = 0;
	char *m_path_buffer = 0;
	float *m_colscores = 0;
	parasail_result_t *m_parasail_result = 0;
	uint8_t *m_parasail_nomalloc_workspace = 0;
	uint m_parasail_nomalloc_workspace_bytes = 0;
	uint *m_pos_is = 0;
	uint *m_pos_js = 0;
	uint *m_considered_vec = 0;
	uint *m_preserved_vec = 0;

private:
	flat_bench2_thread_data() = delete;

public:
	flat_bench2_thread_data(uint nfeat)
		{
		const uint maxL = flat_params::m_maxL;
		m_scratch_rows = myalloc(float, 2*maxL + 3);
		m_scratch_pssms = myalloc(const float *, nfeat);
		m_TB = myalloc(uint8_t, maxL*maxL);
		m_path_buffer = myalloc(char, 2*maxL);
		m_colscores = myalloc(float, 2*maxL);
		m_parasail_result = 0;
		m_parasail_nomalloc_workspace_bytes = (uint)
			parasail_nomalloc_sw_striped_profile_avx2_256_16_workspace_bytes(maxL);
		m_parasail_nomalloc_workspace =
			myalloca(uint8_t, m_parasail_nomalloc_workspace_bytes);
		m_pos_is = myalloc(uint, maxL);
		m_pos_js = myalloc(uint, maxL);
		m_considered_vec = myalloc(uint, maxL);
		m_preserved_vec = myalloc(uint, maxL);
		}

	~flat_bench2_thread_data()
		{
		myfree(m_scratch_pssms);
		myfree(m_scratch_rows);
		myfree(m_colscores);
		myfree(m_pos_is);
		myfree(m_pos_js);
		myfree(m_considered_vec);
		myfree(m_preserved_vec);
		myfree(m_TB);
		myfree(m_path_buffer);
		if (m_parasail_result != 0)
			parasail_result_free(m_parasail_result);
		myfreea(m_parasail_nomalloc_workspace);
		}
	};
