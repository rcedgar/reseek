#pragma once

#include "flat_features.h"

class flat_profiles
	{
public:
	flat_features *m_ff = 0;
	vector<string> m_labels;
	vector<vector<uint8_t> > m_profiles;
	unordered_map<string, uint> m_label2idx;

public:
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
		return m_profiles[i].data();
		}

	uint8_t *get_rev_profile(uint i) const;

	uint32_t get_length(uint i) const
		{
		assert(m_ff != 0);
		assert(m_ff->m_nfeat > 0);
		assert(i < m_profiles.size());
		uint32_t Ln = uint32_t(m_profiles[i].size());
		assert(Ln%m_ff->m_nfeat == 0);
		uint32_t L = Ln/m_ff->m_nfeat;
		return L;
		}

	const string &get_label(uint i) const
		{
		assert(i < m_profiles.size());
		return m_labels[i];
		}

	uint get_nprof() const { return uint(m_profiles.size()); }

	void profile_to_fasta(FILE *f, uint i) const;
	};