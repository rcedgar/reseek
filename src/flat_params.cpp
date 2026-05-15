#include "myutils.h"
#include "flat_params.h"

/////////////////////
// Chain quantization
// MUST RE-TRAIN THRESHOLDS AND LOGODDS
///////////////////////////////////////
uint32_t flat_params::m_nn_min_offset = 12;
uint32_t flat_params::m_distmx_bandwidth = 256;
uint32_t flat_params::m_turnd_w = 5;
uint32_t flat_params::m_angle_n = 4;
////////////////////////////////////

///////////////////////
// LDDT -- special case
float flat_params::m_LDDT_R0 = 15;
static const float thresholds[] = { 0.5, 1, 2, 4 };
const float *flat_params::m_LDDT_thresholds = thresholds;
uint flat_params::m_LDDT_nr_thresholds
	= sizeof(thresholds)/sizeof(thresholds[0]);
///////////////////////////////////////////////

////////////////////////////////////////////////////
// TUNABLE WITHOUT RETRAINING THRESHOLDS AND LOGODDS
// Gap parameters
float flat_params::m_open = FLT_MAX;
float flat_params::m_ext = FLT_MAX;

// Test statistic weights
float flat_params::m_self_w = 0;
float flat_params::m_rev_w = 0;
float flat_params::m_lddt_w = 0;
float flat_params::m_lddtx_w = 0;
float flat_params::m_dali_w = 0;
float flat_params::m_dalix_w = 0;

// Nu filter
float flat_params::m_nu_filter_self_w = 0;
float flat_params::m_nu_filter_rev_w = 0;
float flat_params::m_nu_filter_min_fwd_score = -9999;
float flat_params::m_nu_filter_min_combined_score = -9999;
/////////////////////////////////////////////////////

void flat_params::set_params(
	const vector<string> &names,
	const vector<float> &values)
	{
	assert(names.size() == values.size());
	for (size_t i = 0; i < names.size(); ++i)
		{
		const string &name = names[i];
		float value = values[i];

		if (name == "gap2")
			{
			m_open = value;
			m_ext = value/10;
			}
#define x(param_name, m_name)	else if (name == #param_name) m_name = value;
#include "tunable_flat_params.h"
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
		flat_params::m_lddtx_w > 0;
	}

bool flat_params::need_self()
	{
	return flat_params::m_self_w > 0;
	}

bool flat_params::need_nu_self()
	{
	return flat_params::m_nu_filter_self_w > 0;
	}

bool flat_params::need_reverse()
	{
	return flat_params::m_rev_w > 0;
	}

bool flat_params::need_alignx()
	{
	return
		need_self() ||
		need_reverse();
	}

void flat_params::logme()
	{
	Log("\n");
#define w(x)	Log("%10.3g  %s\n", m_##x, #x)
	w(open);
	w(ext);
	w(self_w);
	w(rev_w);
	w(lddt_w);
	w(lddtx_w);
	w(dali_w);
	w(dalix_w);
	w(nu_filter_self_w);
	w(nu_filter_rev_w);
#undef w

#define w(x)	Log("%10u  %s\n", m_##x, #x)
	w(nn_min_offset);
	w(distmx_bandwidth);
	w(turnd_w);
	w(angle_n);
	w(nu_filter_min_fwd_score);
	w(nu_filter_min_combined_score);
#undef x

	Log("LDDT: R0=%.3g thresholds", m_LDDT_R0);
	for (uint i = 0; i < m_LDDT_nr_thresholds; ++i)
		Log(" %.1f", m_LDDT_thresholds[i]);
	Log("\n");
	}
