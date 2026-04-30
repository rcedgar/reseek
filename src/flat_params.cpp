#include "myutils.h"
#include "flat_params.h"

float flat_params::m_open = FLT_MAX;
float flat_params::m_ext = FLT_MAX;
float flat_params::m_self_w;
float flat_params::m_rev_w;
float flat_params::m_lddt_w;
float flat_params::m_dali_w;
float flat_params::m_dalix_w;
float flat_params::m_entropy_w;

bool flat_params::m_oldts = false;
float flat_params::m_oldts_dpw;
float flat_params::m_oldts_lddtw;
float flat_params::m_oldts_revtsw;
float flat_params::m_oldts_ladd;

uint32_t flat_params::m_nn_min_offset = 12;
uint32_t flat_params::m_distmx_bandwidth = 256;
uint32_t flat_params::m_turnd_w = 5;
uint32_t flat_params::m_angle_n = 4;

void flat_params::set_params(
	const vector<string> &names,
	const vector<float> &values)
	{
	m_oldts = false;
	assert(names.size() == values.size());
	for (size_t i = 0; i < names.size(); ++i)
		{
		const string &name = names[i];
		float value = values[i];
		if (name == "open")
			m_open = value;
		else if (name == "ext")
			m_ext = value;
		else if (name == "selfw")
			m_self_w = value;
		else if (name == "revw")
			m_rev_w = value;
		else if (name == "dali")
			m_dali_w = value;
		else if (name == "dalix")
			m_dalix_w = value;
		else if (name == "lddt")
			m_lddt_w = value;
		else if (name == "entropy")
			m_entropy_w = value;
		else if (name == "gap2")
			{
			m_open = value;
			m_ext = value/10;
			}
		else if (name == "oldts_dpw")
			{
			m_oldts = true;
			m_oldts_dpw = value;
			}
		else if (name == "oldts_revtsw")
			{
			m_oldts = true;
			m_oldts_revtsw = value;
			}
		else if (name == "oldts_lddtw")
			{
			m_oldts = true;
			m_oldts_lddtw = value;
			}
		else if (name == "oldts_ladd")
			{
			m_oldts = true;
			m_oldts_ladd = value;
			}
		else
			Die("flat_params::setparams() %s=%.3g",
				name.c_str(), value);
		}
	asserta(m_open != FLT_MAX);
	asserta(m_ext != FLT_MAX);
	}

bool flat_params::need_distmx()
	{
	return
		flat_params::m_dalix_w > 0 ||
		flat_params::m_dali_w > 0 ||
		flat_params::m_lddt_w > 0 ||
		flat_params::m_entropy_w > 0 ||
		flat_params::m_oldts_lddtw > 0;
	}

bool flat_params::need_prof()
	{
	return flat_params::m_entropy_w > 0;
	}

bool flat_params::need_self()
	{
	return flat_params::m_self_w > 0 || flat_params::m_oldts_revtsw;
	}

bool flat_params::need_reverse()
	{
	return flat_params::m_rev_w > 0;
	}

bool flat_params::need_alignx()
	{
	return
		need_prof() ||
		need_self() ||
		need_reverse();
	}
