#pragma once

#include "hitdata.h"

struct reseek_partial_entry
	{
	reseek_hit hit;
	bool nu_only = false;
	};

struct parent_aln_span
	{
	uint32_t parent_lo_q = 0;
	uint32_t parent_hi_q = 0;
	uint32_t parent_lo_t = 0;
	uint32_t parent_hi_t = 0;
	};

void reseek_hit_parent_aln_span(const reseek_hit &hit, parent_aln_span &span);

float reseek_hit_combined_ts_stub(
	const vector<reseek_partial_entry> &partials,
	uint chosen_index,
	const vector<parent_aln_span> &spans);

reseek_hit reseek_hit_merge_partials(const vector<reseek_partial_entry> &partials);
