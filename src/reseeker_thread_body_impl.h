#pragma once

#include "hitdata.h"

struct nu_pass
	{
	uint qidx = UINT_MAX;
	uint j = UINT_MAX;
	int nu_fwd_score = 0;
	int nu_rev_score = 0;
	float nu_combined_score = 0;
	};

void reseeker_thread_body_impl(uint threadidx, bool use_nusort);
