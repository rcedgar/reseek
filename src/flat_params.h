#pragma once

class flat_params
	{
public:
	static float m_open;
	static float m_ext;
	static float m_self_w;
	static float m_rev_w;
	static float m_lddt_w;
	static float m_dali_w;
	static float m_dalix_w;
	static float m_entropy_w;

	static bool m_oldts;
	static float m_oldts_dpw;
	static float m_oldts_lddtw;
	static float m_oldts_revtsw;
	static float m_oldts_ladd;

	static uint32_t m_nn_w;

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