#pragma once

class tdata
	{
public:
	static tdata *m_tdatavec;
	static uint32_t m_nt;

public:
	bool m_isq = false;
	const string m_label;
	const flat_chain_t *m_chain = 0;
	uint m_L = 0;
	sid_t *m_distmx = 0;
	uint8_t *m_codeseq_nu = 0;
	uint8_t *m_mega_prof = 0;

public:
	tdata()
		{
		m_isq = false;
		}
	};
