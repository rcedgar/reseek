#pragma once

class flat_chain_t;

class hitdata
	{
public:
	const flat_chain_t *query = 0;
	const flat_chain_t *target = 0;
	const char *path = 0;
	uint ncol = 0;
	float nu_self_score = 0;
	float mega_self_score = 0;
	float nu_fwd_score = 0;
	float nu_rev_score = 0;
	float mega_fwd_score = 0;
	float mega_rev_score = 0;

	float lddt = 0;
	float dali = 0;

	float TS_fold = 0;
	float TS_sf = 0;
	float TS_fam = 0;

public:
	void reset()
		{
		query = 0;
		target = 0;
		path = 0;
		ncol = 0;
		nu_self_score = 0;
		mega_self_score = 0;
		nu_fwd_score = 0;
		nu_rev_score = 0;
		mega_fwd_score = 0;
		mega_rev_score = 0;
		lddt = 0;
		dali = 0;
		TS_fold = 0;
		TS_sf = 0;
		TS_fam = 0;
		}
	};
