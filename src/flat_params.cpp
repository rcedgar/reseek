#include "myutils.h"
#include "flat_params.h"

// Gap parameters
float flat_params::m_open = FLT_MAX;
float flat_params::m_ext = FLT_MAX;

// Test statistic weights
float flat_params::m_self_w;
float flat_params::m_rev_w;
float flat_params::m_lddt_w;
float flat_params::m_lddtx_w;
float flat_params::m_dali_w;
float flat_params::m_dalix_w;

/***
Nu self & rev weighting
=======================
Optimal parameters for compound aa4+pm2+sec32, fixed weights from above, fwd&rev:
	m_Scores[HitIdx] = Score_fwd - RevWeight*Score_rev -
		SelfWeight*(SelfScore_rev_i + SelfScore_rev_j);
	numegarev.log:05:42 527Mb  >>>0.00016[1.38512] latinclimb:HJ2/2:explore+selfw /1.30/ selfw=5.0E-01;revw=2.7E-01;
	=>selfw=0.5;revw=0.27;
***/
// Nu filter
// selfw=5.0E-01;revw=2.7E-01;
float flat_params::m_nu_filter_self_w = 0.5f;
float flat_params::m_nu_filter_rev_w = 0.27f;
int flat_params::m_min_nu_fwd_score = 130;
int flat_params::m_min_nu_combined_score = 40;

// Chain quantization
uint32_t flat_params::m_nn_min_offset = 12;
uint32_t flat_params::m_distmx_bandwidth = 256;
uint32_t flat_params::m_turnd_w = 5;
uint32_t flat_params::m_angle_n = 4;

float flat_params::m_LDDT_R0 = 15;
static const float thresholds[] = { 0.5, 1, 2, 4 };
const float *flat_params::m_LDDT_thresholds = thresholds;
uint flat_params::m_LDDT_nr_thresholds
	= sizeof(thresholds)/sizeof(thresholds[0]);

void flat_params::set_params(
	const vector<string> &names,
	const vector<float> &values)
	{
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
		else if (name == "lddtx")
			m_lddtx_w = value;
		else if (name == "nfselfw")
			m_nu_filter_self_w = value;
		else if (name == "nfrevw")
			m_nu_filter_rev_w = value;
		else if (name == "gap2")
			{
			m_open = value;
			m_ext = value/10;
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
	w(min_nu_fwd_score);
#undef x

	Log("LDDT: R0=%.3g thresholds", m_LDDT_R0);
	for (uint i = 0; i < m_LDDT_nr_thresholds; ++i)
		Log(" %.1f", m_LDDT_thresholds[i]);
	Log("\n");
	}
