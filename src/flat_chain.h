#pragma once

#include "flat_base.h"

static const size_t RESERVE_CHAIN_LENGTH = 400;

class flat_chain
	{
public:
	string m_label;
	chainxyz_t *m_xyz;
	chainaa_t *m_aa;
	vector<string> m_lines;

	flat_chain()
		{
		m_label = "(undef)";
		m_xyz = nullptr;
		m_aa = nullptr;
		}

	flat_chain(
		const string &label,
		const vector<char> aas,
		const vector<float> &Xs,
		const vector<float> &Ys,
		const vector<float> &Zs);

	~flat_chain()
		{
		clear();
		}

	void clear()
		{
		m_lines.clear();
		release0(m_xyz);
		release0(m_aa);
		}

	void set_xyz(const vector<float> &Xs,
		const vector<float> &Ys, const vector<float> &Zs);
	void set_aa(const vector<char> &aas);
	bool from_pdb_lines(const string &label,
		const vector<string> &lines, bool save_lines);
	void to_fasta(const string &fn) const;
	void to_fasta(FILE *f) const;

	uint32_t get_length() const { return m_aa == 0 ? 0 : m_aa->m_size; }
	};