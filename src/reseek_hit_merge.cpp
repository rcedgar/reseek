#include "myutils.h"
#include "reseek_hit_merge.h"
#include "cigar.h"

void reseek_hit_parent_aln_span(const reseek_hit &hit, parent_aln_span &span)
	{
	uint qL = 0;
	uint tL = 0;
	if (!hit.cigar.empty())
		CIGARToLs(hit.cigar, qL, tL);
	span.parent_lo_q = hit.query.parent_lo + hit.fwd.lo_q;
	span.parent_lo_t = hit.target.parent_lo + hit.fwd.lo_t;
	span.parent_hi_q = span.parent_lo_q + qL;
	span.parent_hi_t = span.parent_lo_t + tL;
	}

static uint best_partial_index(const vector<reseek_partial_entry> &partials)
	{
	asserta(!partials.empty());
	uint best_i = 0;
	for (uint i = 1; i < uint(partials.size()); ++i)
		{
		const reseek_partial_entry &a = partials[best_i];
		const reseek_partial_entry &b = partials[i];
		float score_a = a.nu_only ? a.hit.nu_combined_score : a.hit.TS;
		float score_b = b.nu_only ? b.hit.nu_combined_score : b.hit.TS;
		if (score_b > score_a)
			best_i = i;
		}
	return best_i;
	}

float reseek_hit_combined_ts_stub(
	const vector<reseek_partial_entry> &partials,
	uint chosen_index,
	const vector<parent_aln_span> &spans)
	{
	(void) spans;
	asserta(chosen_index < partials.size());
	const reseek_partial_entry &chosen = partials[chosen_index];
	if (chosen.nu_only)
		return chosen.hit.nu_combined_score;
	return chosen.hit.TS;
	}

reseek_hit reseek_hit_merge_partials(const vector<reseek_partial_entry> &partials)
	{
	asserta(!partials.empty());
	const uint chosen_index = best_partial_index(partials);
	const reseek_partial_entry &chosen = partials[chosen_index];

	vector<parent_aln_span> spans(partials.size());
	for (uint i = 0; i < uint(partials.size()); ++i)
		reseek_hit_parent_aln_span(partials[i].hit, spans[i]);

	reseek_hit merged = chosen.hit;
	merged.TS = reseek_hit_combined_ts_stub(partials, chosen_index, spans);
	return merged;
	}
