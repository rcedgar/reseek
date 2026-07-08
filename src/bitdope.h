#pragma once

#include "lookup.h"

static const uint32_t MAGIC	= 0xd05e;

class bitdope
	{
public:
	uint8_t *m_dope = 0;
	lookup *m_look = 0;
	uint m_nhit = 0;
	uint m_ndom = 0;
	vector<bool> m_square;

public:
	void from_file(const string &fn)
		{
		//asserta(m_look != 0);

		m_nhit = 0;
		uint32_t magic;
		FILE *f = OpenStdioFile(fn);
		ReadStdioFile(f, &magic, sizeof(magic));
		asserta(magic == MAGIC);
		ReadStdioFile(f, &m_ndom, sizeof(m_ndom));
		//asserta(m_ndom == m_look->get_ndom());
		uint32_t K = triangle_get_K(m_ndom);
		uint32_t bytes = (K + 7)/8;
		m_dope = myalloc(uint8_t, bytes);

		ReadStdioFile(f, m_dope, bytes);
		ReadStdioFile(f, &magic, sizeof(magic));
		asserta(magic == MAGIC);
		CloseStdioFile(f);

		for (uint i = 0; i < bytes; ++i)
			{
			uint8_t b = m_dope[i];
			for (uint j = 0; j < 8; ++j)
				{
				if (b & (1 << j))
					++m_nhit;
				}
			}
		}

	bool in_dope_k(uint k) const
		{
		if (m_dope == 0) return true;
		byte b = m_dope[k/8];
		return b & (1 << k%8);
		}
	bool in_dope_ij(uint i, uint j) const
		{
		if (m_dope == 0) return true;
		if (i == j) return false;
		uint k = triangle_ij_to_k2(i, j, m_ndom);
		return in_dope_k(k);
		}

	void set_square()
		{
		m_square.clear();
		m_square.resize(m_ndom*m_ndom);
		uint n = 0;
		for (uint i = 0; i < m_ndom; ++i)
			{
			for (uint j = 0; j < i; ++j)
				{
				if (in_dope_ij(i, j))
					{
					m_square[i*m_ndom + j] = true;
					m_square[j*m_ndom + i] = true;
					n += 2;
					}
				}
			}
		if (n != 2*m_nhit)
			Die("n=%u, 2*m_nhit=%u", n, 2*m_nhit);
		}

	bool in_square_ij(uint i, uint j) const
		{
		assert(i < m_ndom);
		assert(j < m_ndom);
		return m_square[i*m_ndom + j];
		}
	};