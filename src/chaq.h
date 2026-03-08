#pragma once

#include "flat_base.h"
#include "flat_chain.h"

/***
Chain quantizer
Data derived from a chain is created here.
The chain itself *m_chain (label+aa+xyz) owned elsewhere.
Caller is responsible for ensuring that m_chain remains valid.
Data is created on demand and persists until explicitly freed
	or chaq::init() is called for a second time.
Pointers are private to avoid trivial copying.
External caller must ensure that down(p) is called.
Member functions of chaq access member pointers directly, rarely /
	never calls down() on m_* ptr to ensure that data is cached.
***/
class chaq
	{
public:
	const flat_chain *m_chain = 0;

private:
	////////////////////////////
	// Cached data
	nnvec_t *m_pen = 0;
	nnvec_t *m_men = 0;
	floatvec_t *m_pendist = 0;
	floatvec_t *m_mendist = 0;
	ss3_t *m_ss3 = 0;
	megaprof_t *m_megaprof = 0;
	chaindistmx_t *m_distmx = 0;
	////////////////////////////

public:
	void clear()
		{
		down0(m_pen);
		down0(m_men);
		down0(m_ss3);
		down0(m_megaprof);
		down0(m_distmx);
		}

	void init(const flat_chain *chain)
		{
		clear();
		m_chain = chain;
		}

	uint32_t get_length() const { assert(m_chain); return m_chain->get_length(); }
	const nnvec_t *get_pen();
	const nnvec_t *get_men();
	const ss3_t *get_ss3();
	const megaprof_t *get_megaprof();
	const chaindistmx_t *get_distmx(uint M);
	uint16_t get_nen(uint i) const;
	uint16_t get_ren(uint i) const;
	uint16_t get_pen(uint i) const;
	uint16_t get_men(uint i) const;
	float get_pen_dist_float(uint i) const;
	uint16_t get_pen_dist_ic(uint i) const;
	};
