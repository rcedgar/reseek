#pragma once

#include "triangle.h"

enum LOOK_TRUTH
	{
	LT_Undef,
	LT_SAME_FAM,
	LT_SAME_SF,
	LT_DIFF_FAM_SAME_SF,
	LT_DIFF_SF_SAME_FOLD,
	LT_SAME_FOLD
	};

class flat_chain_t;

class lookup
	{
public:
	vector<string> m_doms;
	vector<string> m_fams;
	vector<string> m_sfs;
	vector<string> m_folds;
	unordered_map<string, uint> m_dom2idx;
	unordered_map<string, uint> m_fam2idx;
	unordered_map<string, uint> m_sf2idx;
	unordered_map<string, uint> m_fold2idx;
	vector<uint> m_domidx2famidx;
	vector<uint> m_domidx2sfidx;
	vector<uint> m_domidx2foldidx;
	vector<uint> m_famidx2ndom;
	vector<uint> m_sfidx2ndom;
	vector<uint> m_foldidx2ndom;
	uint m_NT = 0;
	uint m_NF = 0;
	uint m_NI = 0;
	uint m_pair_count = 0;
	LOOK_TRUTH m_LT = LT_SAME_SF;

public:
	lookup()
		{
		m_LT = LT_SAME_SF;
		if (optset_truth)
			{
			const string &t = opt(truth);
			if (t == "fam")
				m_LT = LT_SAME_FAM;
			else if (t == "sf")
				m_LT = LT_SAME_SF;
			else if (t == "sfx")
				m_LT = LT_DIFF_FAM_SAME_SF;
			else if (t == "fold")
				m_LT = LT_SAME_FOLD;
			else if (t == "foldx")
				m_LT = LT_DIFF_SF_SAME_FOLD;
			else
				Die("Invalid -truth '%s'", t.c_str());
			}
		}

	~lookup()
		{
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
		}

	const char *get_truthstr() const
		{
		if (m_LT == LT_SAME_FAM)
			return "fam";
		else if (m_LT == LT_SAME_SF)
			return "sf";
		else if (m_LT == LT_DIFF_FAM_SAME_SF)
			return "sfx";
		else if (m_LT == LT_SAME_FOLD)
			return "fold";
		else if (m_LT == LT_DIFF_SF_SAME_FOLD)
			return "foldx";
		else
			Die("get_truthstr()");
		return "ERROR";
		}

	uint get_ndom() const { return uint(m_doms.size()); }

	void from_tsv(const string &fn);
	void from_labels(const vector<string> &labels);
	void to_tsv(const string &fn);
	void fill();
	void fill_fam();
	void fill_sfx();
	void fill_sf();
	void fill_fold();
	void fill_foldx();
	void stats();

	uint get_singleton_count() const;

	const string &get_dom(uint domidx) const
		{
		assert(domidx < m_doms.size());
		return m_doms[domidx];
		}

	const char *get_fam(uint domidx) const
		{
		assert(domidx < m_doms.size());
		assert(domidx < m_domidx2famidx.size());
		uint famidx = m_domidx2famidx[domidx];
		assert(famidx < m_fams.size());
		return m_fams[famidx].c_str();
		}
	
	void get_dom_scopid(uint domidx, string &label) const
		{
		assert(domidx < m_doms.size());
		assert(domidx < m_domidx2sfidx.size());
		uint sfidx = m_domidx2sfidx[domidx];
		assert(sfidx < m_sfs.size());
		label = m_doms[domidx] + "/" + m_sfs[sfidx];
		}

	uint get_domidx(const string &lab, bool failok = false) const
		{
		void trunc_label(const string &Label, string &TruncatedLabel);

		string dom;
		trunc_label(lab, dom);
		unordered_map<string, uint>::const_iterator iter =
			m_dom2idx.find(dom);
		if (iter == m_dom2idx.end())
			{
			if (failok)
				return UINT_MAX;
			Die("get_domidx(%s)", lab.c_str());
			}
		return iter->second;
		}

	uint get_pair_count_upper_triangle_with_diagonal() const
		{
		uint K = triangle_get_K(get_ndom());
		return K;
		}

	uint get_pair_idx_upper_triangle_with_diagonal(
		uint domidx1, uint domidx2) const
		{
		uint minidx = min(domidx1, domidx2);
		uint maxidx = max(domidx1, domidx2);
		uint k = triangle_ij_to_k(minidx, maxidx,
			uint(m_doms.size()));
		return k;
		}

	bool same_sf_ij(uint i, uint j) const
		{
		uint sfidx_i = m_domidx2sfidx[i];
		uint sfidx_j = m_domidx2sfidx[j];
		return sfidx_i == sfidx_j;
		}

	bool same_fam_ij(uint i, uint j) const
		{
		uint famidx_i = m_domidx2famidx[i];
		uint famidx_j = m_domidx2famidx[j];
		return famidx_i == famidx_j;
		}

	bool same_sf_k(uint k) const
		{
		uint i, j;
		triangle_k_to_ij(k, uint(m_doms.size()), i, j);
		return same_sf_ij(i, j);
		}

	bool same_fold_ij(uint i, uint j) const
		{
		uint foldidx_i = m_domidx2foldidx[i];
		uint foldidx_j = m_domidx2foldidx[j];
		return foldidx_i == foldidx_j;
		}

	bool same_fold_k(uint k) const
		{
		uint i, j;
		triangle_k_to_ij(k, uint(m_doms.size()), i, j);
		return same_fold_ij(i, j);
		}

	bool is_ignored_ij(uint i, uint j) const
		{
		if (m_LT == LT_DIFF_SF_SAME_FOLD)
			{
			uint foldidx_i = m_domidx2foldidx[i];
			uint foldidx_j = m_domidx2foldidx[j];
			uint sfidx_i = m_domidx2sfidx[i];
			uint sfidx_j = m_domidx2sfidx[j];
			return sfidx_i == sfidx_j && foldidx_i == foldidx_j;
			}
		else if (m_LT == LT_DIFF_FAM_SAME_SF)
			{
			uint famidx_i = m_domidx2famidx[i];
			uint famidx_j = m_domidx2famidx[j];
			return famidx_i == famidx_j;
			}
		else
			return false;
		}

	bool is_ignored_k(uint k) const
		{
		uint i, j;
		triangle_k_to_ij(k, uint(m_doms.size()), i, j);
		return is_ignored_ij(i, j);
		}

	bool is_tp_ij(uint i, uint j) const
		{
		if (m_LT == LT_SAME_SF)
			return same_sf_ij(i, j);
		else if (m_LT == LT_SAME_FAM)
			return same_fam_ij(i, j);
		else if (m_LT == LT_SAME_FOLD)
			return same_fold_ij(i, j);
		else if (m_LT == LT_DIFF_FAM_SAME_SF)
			return !same_fam_ij(i, j) && same_sf_ij(i, j);
		else if (m_LT == LT_DIFF_SF_SAME_FOLD)
			return !same_sf_ij(i, j) && same_fold_ij(i, j);
		return false;
		}

	bool is_tp_k(uint k) const
		{
		uint i, j;
		triangle_k_to_ij(k, uint(m_doms.size()), i, j);
		return is_tp_ij(i, j);
		}

	void sort_chains(
		const vector<flat_chain_t *> &chains,
		vector<flat_chain_t *> &sorted_chains) const;
	};
