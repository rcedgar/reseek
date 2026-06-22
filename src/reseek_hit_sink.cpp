#include "myutils.h"
#include "reseek_hit_sink.h"
#include "reseek_hit_merge.h"
#include "reseeker.h"
#include <map>

static vector<vector<reseek_partial_entry> > *s_thread_buffers = 0;
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

void reseek_hit_sink_begin(uint thread_count)
	{
	asserta(thread_count > 0);
	if (s_thread_buffers == 0)
		s_thread_buffers = new vector<vector<reseek_partial_entry> >();
	s_thread_buffers->clear();
	s_thread_buffers->resize(thread_count);
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

void reseek_hit_sink_thread_end(uint threadidx)
	{
	if (s_thread_partial == 0 || s_thread_partial->empty())
		return;
	asserta(s_thread_buffers != 0);
	asserta(threadidx < s_thread_buffers->size());
	vector<reseek_partial_entry> &dest = (*s_thread_buffers)[threadidx];
	const uint n = uint(s_thread_partial->size());
	for (uint i = 0; i < n; ++i)
		dest.push_back((*s_thread_partial)[i]);
	s_thread_partial->clear();
	}

static void collect_all_partials(vector<reseek_partial_entry> &all)
	{
	all.clear();
	if (s_thread_buffers == 0)
		return;
	const uint nthread = uint(s_thread_buffers->size());
	for (uint t = 0; t < nthread; ++t)
		{
		const vector<reseek_partial_entry> &tb = (*s_thread_buffers)[t];
		const uint n = uint(tb.size());
		for (uint i = 0; i < n; ++i)
			all.push_back(tb[i]);
		}
	for (uint t = 0; t < nthread; ++t)
		(*s_thread_buffers)[t].clear();
	}

void reseek_hit_sink_flush(FILE *fhit)
	{
	vector<reseek_partial_entry> all;
	collect_all_partials(all);
	if (all.empty())
		return;

	map<parent_pair_key, vector<reseek_partial_entry> > groups;
	const uint n = uint(all.size());
	for (uint i = 0; i < n; ++i)
		{
		const reseek_partial_entry &e = all[i];
		parent_pair_key k = reseek_hit_parent_pair_key(e.hit);
		groups[k].push_back(e);
		}

	for (map<parent_pair_key, vector<reseek_partial_entry> >::const_iterator
			iter = groups.begin(); iter != groups.end(); ++iter)
		{
		const vector<reseek_partial_entry> &partials = iter->second;
		reseek_hit merged = reseek_hit_merge_partials(partials);
		const bool nu_only = partials[0].nu_only;
		reseek_hit_emit_tsv(fhit, nu_only, merged);
		}
	}
