#pragma once

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
