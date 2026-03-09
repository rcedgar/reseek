#include "myutils.h"
#include "flat_base.h"
#include "chaq.h"
#include "flat_distmx.h"

const chaindistmx_t *chaq::get_distmx(uint M)
	{
	assert(m_chain);
	if (m_distmx) return m_distmx;
	const uint32_t L = get_length();
	_chkmem();//@@
	m_distmx = create_chaindistmx(L, M);
	_chkmem();//@@
	const uint16_t *xyz = m_chain->m_xyz->m_data;
	uint16_t *distmx = m_distmx->m_data;
	fill_flat_distmx(xyz, L, M, distmx);
	up(m_distmx);
	return m_distmx;
	}

float chaq::get_pen_dist_float(uint i) const
	{
	return 0;//@@TODO
	}

uint16_t chaq::get_pen_dist_ic(uint i) const
	{
	uint16_t nen = get_pen(i);
	uint16_t *distmx = m_distmx->m_data;
	//uint32_t k = band_ij_to_k(i, nen);
	return 0;//@@TODO
	}

uint16_t chaq::get_nen(uint i) const
	{
	assert(false);
	return 0;
	}

uint16_t chaq::get_pen(uint i) const
	{
	assert(m_pen);
	assert(i < m_pen->m_size);
	return m_pen->m_data[i];
	}

uint16_t chaq::get_men(uint i) const
	{
	assert(m_men);
	assert(i < m_men->m_size);
	return m_men->m_data[i];
	}
