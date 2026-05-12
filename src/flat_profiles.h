#pragma once

#include "flat_alphas.h"
#include "lookup.h"

class flat_profiles
	{
public:
	vector<string> m_labels;
	vector<uint> m_lengths;
	vector<uint8_t *> m_profiles;
	unordered_map<string, uint> m_label2idx;
	vector<uint8_t *> m_nu_codeseqs;

public:
	void from_chains(const vector<flat_chain_t *> &chains);
	void from_chains_lookup(const lookup &look,
		const vector<flat_chain_t *> &chains);

	void write_nu_hexfasta(const string &fn) const;

	void read_profiles_faprof(
		const string &faproffn,
		vector<string> &feature_names);

	void read_profiles_from_fastas(
		const vector<string> &fafns,
		const unordered_map<string, uint> &label2idx);

	void check_profiles() const;
	void check_profile(uint i) const;

	const uint8_t *get_profile(uint i) const
		{
		assert(i < m_profiles.size());
		return m_profiles[i];
		}

	uint8_t *get_rev_profile(uint i) const;

	uint32_t get_length(uint i) const
		{
		assert(flat_alphas::m_nfeat > 0);
		assert(i < m_profiles.size());
		uint32_t L = uint32_t(m_lengths[i]);
		return L;
		}

	const string &get_label(uint i) const
		{
		assert(i < m_profiles.size());
		return m_labels[i];
		}

	uint get_nprof() const { return uint(m_profiles.size()); }

	void profile_to_fasta(FILE *f, uint i) const;

	uint8_t *make_profile(const flat_chain_t &chain) const;

	void set_nu_codeseqs(const string &hexfastafn = "");
	void set_nu_codeseq(uint fi_aa20, uint fi_pm2, uint fi_sec32, uint idx);
	};