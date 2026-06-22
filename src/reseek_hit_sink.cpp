#include "myutils.h"
#include "reseek_hit_sink.h"
#include "reseeker.h"
#include <map>

struct reseek_partial_entry
	{
	reseek_hit hit;
	bool nu_only = false;
	};

static vector<reseek_partial_entry> *s_buffered = 0;
static thread_local vector<reseek_partial_entry> *s_thread_partial = 0;

static vector<reseek_partial_entry> &thread_partials()
	{
	if (s_thread_partial == 0)
		s_thread_partial = new vector<reseek_partial_entry>();
	return *s_thread_partial;
	}

static bool hit_is_sliced(const reseek_hit &hit)
	{
	return hit.query.is_sliced() || hit.target.is_sliced();
	}

void reseek_hit_sink_begin()
	{
	if (s_buffered == 0)
		s_buffered = new vector<reseek_partial_entry>();
	s_buffered->clear();
	}

void reseek_hit_sink_submit(const reseek_hit &hit, bool nu_only)
	{
	if (!hit_is_sliced(hit))
		{
		reseek_hit_emit_tsv(reseeker::m_fhit, nu_only, hit);
		return;
		}

	reseek_partial_entry e;
	e.hit = hit;
	e.nu_only = nu_only;
	thread_partials().push_back(e);
	}

void reseek_hit_sink_thread_end()
	{
	if (s_thread_partial == 0 || s_thread_partial->empty())
		return;
	asserta(s_buffered != 0);
	const uint n = uint(s_thread_partial->size());
	for (uint i = 0; i < n; ++i)
		s_buffered->push_back((*s_thread_partial)[i]);
	s_thread_partial->clear();
	}

static const reseek_partial_entry &best_partial_entry(
	const vector<reseek_partial_entry> &entries)
	{
	asserta(!entries.empty());
	uint best_i = 0;
	for (uint i = 1; i < uint(entries.size()); ++i)
		{
		const reseek_partial_entry &a = entries[best_i];
		const reseek_partial_entry &b = entries[i];
		float score_a = a.nu_only ? a.hit.nu_combined_score : a.hit.TS;
		float score_b = b.nu_only ? b.hit.nu_combined_score : b.hit.TS;
		if (score_b > score_a)
			best_i = i;
		}
	return entries[best_i];
	}

void reseek_hit_sink_flush(FILE *fhit)
	{
	if (s_buffered == 0 || s_buffered->empty())
		return;

	map<parent_pair_key, vector<reseek_partial_entry> > groups;
	const uint n = uint(s_buffered->size());
	for (uint i = 0; i < n; ++i)
		{
		const reseek_partial_entry &e = (*s_buffered)[i];
		parent_pair_key k = reseek_hit_parent_pair_key(e.hit);
		groups[k].push_back(e);
		}

	for (map<parent_pair_key, vector<reseek_partial_entry> >::const_iterator
			iter = groups.begin(); iter != groups.end(); ++iter)
		{
		const reseek_partial_entry &best = best_partial_entry(iter->second);
		reseek_hit_emit_tsv(fhit, best.nu_only, best.hit);
		}

	s_buffered->clear();
	}
