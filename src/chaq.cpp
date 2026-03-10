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

// 0=helix 1=strand 2=turn 3=loop
// Method from sec_str() in TMalign.cpp Zhang & Skolnick 2005
uint8_t chaq::get_ss4(const sid_t *distmx, uint M, uint L, uint pos) const
	{
	if (pos < 2 || pos + 2 >= L)
		return 3;

	float dis13 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos-2, pos)]);
	float dis14 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos-2, pos+1)]);
	float dis15 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos-2, pos+2)]);
	float dis24 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos-1, pos+1)]);
	float dis25 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos-1, pos+2)]);
	float dis35 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos, pos+2)]);

	const float DH = 2.1;
	if (fabs(dis15 - 6.37) < DH && fabs(dis14 - 5.18) < DH &&
		fabs(dis25 - 5.18) < DH && fabs(dis13 - 5.45) < DH &&
		fabs(dis24 - 5.45) < DH && fabs(dis35 - 5.45) < DH)
		return 0;	// helix

	const float DS = 1.42;
	if (fabs(dis15 - 13) < DS && fabs(dis14 - 10.4) < DS &&
		fabs(dis25 - 10.4) < DS && fabs(dis13 - 6.1) < DS &&
		fabs(dis24 - 6.1) < DS && fabs(dis35 - 6.1) < DS)
		return 1;	// strand

	if (dis15 < 8.2)
		return 2;	// turn

	return 3; // loop
	}

// 0=helix 1=strand 2=other
uint8_t chaq::get_ss3(const sid_t *distmx, uint M, uint L, uint pos) const
	{
	uint8_t ss4 = get_ss4(distmx, M, L, pos);
	return ss4 <= 2 ? ss4 : 2;
	}

void chaq::get_ss4_str(const sid_t *distmx, uint M, uint L, string &ss) const
	{
	ss.clear();
	ss.reserve(L);
	for (uint pos = 0; pos < L; ++pos)
		{
		uint8_t letter = get_ss4(distmx, M, L, pos);
		assert(letter < 4);
		ss += "hst~"[letter];
		}
	}

void chaq::get_ss4_intseq(const sid_t *distmx, uint M, uint L, uint8_t *intseq) const
	{
	for (uint pos = 0; pos < L; ++pos)
		intseq[pos] = get_ss4(distmx, M, L, pos);
	}
