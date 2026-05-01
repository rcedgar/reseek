#pragma once

class flat_params
	{
public:
	static float m_open;
	static float m_ext;
	static float m_self_w;
	static float m_rev_w;
	static float m_lddt_w;
	static float m_lddtx_w;
	static float m_lddtpow_w;
	static float m_dali_w;
	static float m_dalix_w;
	static float m_entropy_w;
	static float m_rotfreetm_w;

	static uint32_t m_nn_min_offset;
	static uint32_t m_distmx_bandwidth;
	static uint32_t m_turnd_w;
	static uint32_t m_angle_n;

	static float m_LDDT_R0;
	static const float *m_LDDT_thresholds;
	static uint m_LDDT_nr_thresholds;

public:
	static void set_params(
		const vector<string> &names,
		const vector<float> &values);

	static bool need_reverse();
	static bool need_distmx();
	static bool need_prof();
	static bool need_self();
	static bool need_alignx();
	};