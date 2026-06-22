#pragma once

#include "myutils.h"

class flat_chain_t;

struct chain_slice
	{
	uint32_t storage_idx = UINT32_MAX;
	uint32_t parent_storage_idx = UINT32_MAX;
	uint32_t parent_lo = 0;
	uint32_t parent_len = 0;
	uint16_t win_id = 0;
	uint16_t slice_len = 0;

	bool is_sliced() const
		{
		if (win_id != 0)
			return true;
		if (parent_lo != 0)
			return true;
		if (parent_len > 0 && slice_len < parent_len)
			return true;
		return false;
		}
	};

struct parent_pair_key
	{
	uint32_t parent_q = UINT32_MAX;
	uint32_t parent_t = UINT32_MAX;

	bool operator<(const parent_pair_key &rhs) const
		{
		if (parent_q != rhs.parent_q)
			return parent_q < rhs.parent_q;
		return parent_t < rhs.parent_t;
		}
	};

struct aln_segment
	{
	uint32_t lo_q = 0;
	uint32_t lo_t = 0;
	uint32_t ncol = 0;
	string path;
	};

struct reseek_hit
	{
	chain_slice query;
	chain_slice target;

	string parent_label_q;
	string parent_label_t;

	aln_segment fwd;
	string cigar;
	string seq_q;
	string seq_t;

	float nu_fwd_score = 0;
	float nu_rev_score = 0;
	float nu_combined_score = 0;
	float mega_fwd_score = 0;
	float mega_rev_score = 0;
	float lddt = 0;
	float dali = 0;
	float TS = 0;
	uint kappa_diag_score = 0;

public:
	void reset()
		{
		query = chain_slice();
		target = chain_slice();
		parent_label_q.clear();
		parent_label_t.clear();
		fwd = aln_segment();
		cigar.clear();
		seq_q.clear();
		seq_t.clear();
		nu_fwd_score = 0;
		nu_rev_score = 0;
		nu_combined_score = 0;
		mega_fwd_score = 0;
		mega_rev_score = 0;
		lddt = 0;
		dali = 0;
		TS = 0;
		kappa_diag_score = 0;
		}

	void set_fwd_path(const char *path, uint ncol, uint lo_q, uint lo_t)
		{
		fwd.lo_q = lo_q;
		fwd.lo_t = lo_t;
		fwd.ncol = ncol;
		if (path == 0 || ncol == 0)
			fwd.path.clear();
		else
			fwd.path.assign(path, ncol);
		}
	};

inline chain_slice chain_slice_identity(uint storage_idx, uint L)
	{
	chain_slice s;
	s.storage_idx = storage_idx;
	s.parent_storage_idx = storage_idx;
	s.parent_lo = 0;
	s.parent_len = L;
	s.win_id = 0;
	s.slice_len = uint16_t(L);
	return s;
	}

inline parent_pair_key reseek_hit_parent_pair_key(const reseek_hit &hit)
	{
	parent_pair_key k;
	k.parent_q = hit.query.parent_storage_idx;
	k.parent_t = hit.target.parent_storage_idx;
	return k;
	}

string flat_chain_aa_string(const flat_chain_t *chain, uint L);
void reseek_hit_fill_aa_cigar(
	reseek_hit &hit,
	const flat_chain_t *query_chain,
	const flat_chain_t *target_chain);

void reseek_hit_emit_tsv(FILE *f, bool nu_only, const reseek_hit &hit);
void reseek_hit_emit_aln_aa(FILE *f, const reseek_hit &hit);
