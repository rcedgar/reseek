#pragma once

#include "triangle.h"

class lookup
	{
public:
	vector<string> m_doms;
	vector<string> m_sfs;
	unordered_map<string, uint> m_dom2idx;
	unordered_map<string, uint> m_sf2idx;
	vector<uint> m_domidx2sfidx;
	vector<uint> m_sfidx2ndom;
	bool *m_tpvec = 0;
	uint m_NT = 0;
	uint m_NF = 0;
	uint m_pair_count = 0;

public:
	~lookup()
		{
		myfree(m_tpvec);
		}

	void reserve()
		{
		const size_t NDOM = 1200;
		const size_t NSF = 1000;

		m_doms.reserve(NDOM);
		m_sfs.reserve(NSF);
		m_dom2idx.reserve(NDOM);
		m_sf2idx.reserve(NSF);
		m_domidx2sfidx.reserve(NDOM);
		m_sfidx2ndom.reserve(NSF);
		}

	void clear()
		{
		m_doms.clear();
		m_sfs.clear();
		m_dom2idx.clear();
		m_sf2idx.clear();
		m_domidx2sfidx.clear();
		m_sfidx2ndom.clear();
		myfree(m_tpvec);
		}

	uint get_ndom() { return uint(m_doms.size()); }

	void from_tsv(const string &fn);
	void from_labels(const vector<string> &labels);
	void to_tsv(const string &fn);
	void fill();
	void stats();

	uint get_pair_count_upper_triangle_with_diagonal()
		{
		uint K = triangle_get_K(get_ndom());
		return K;
		}

	uint get_pair_idx_upper_triangle_with_diagonal(
		uint domidx1, uint domidx2)
		{
		uint minidx = min(domidx1, domidx2);
		uint maxidx = max(domidx1, domidx2);
		uint k = triangle_ij_to_k(minidx, maxidx,
			uint(m_doms.size()));
		return k;
		}

	bool is_tp_ij(uint domidx1, uint domidx2)
		{
		assert(domidx1 < m_domidx2sfidx.size());
		assert(domidx2 < m_domidx2sfidx.size());
		uint sfidx1 = m_domidx2sfidx[domidx1];
		uint sfidx2 = m_domidx2sfidx[domidx2];
		return sfidx1 == sfidx2;
		}

	bool is_tp_k(uint k)
		{
		assert(k < m_pair_count);
		return m_tpvec[k];
		}
	};
