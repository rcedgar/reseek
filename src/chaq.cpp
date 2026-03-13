#include "myutils.h"
#include "flat_base.h"
#include "chaq.h"
#include "flat_distmx.h"

void chaq::create_distmx(const flat_chain_t &chain, 
	chaindistmx_t*& dm, uint M)
	{
	const uint32_t L = chain.get_length();
	dm = chaindistmx_t::newflat(L, M);
	uint16_t *distmx = dm->m_data;
	fill_flat_distmx(chain.m_xyz->m_data, L, M, distmx);
	}

// 0=helix 1=strand 2=turn 3=loop
// Method from sec_str() in TMalign.cpp Zhang & Skolnick 2005
uint8_t chaq::get_ss4(const sid_t *distmx, uint M, uint L, uint pos)
	{
	if (pos < 2 || pos + 2 >= L)
		return 3;

	float dis13 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos-2, pos)]);
	float dis14 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos-2, pos+1)]);
	float dis15 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos-2, pos+2)]);
	float dis24 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos-1, pos+1)]);
	float dis25 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos-1, pos+2)]);
	float dis35 = sid2dist(distmx[banded_i_lt_j_to_k(M, pos, pos+2)]);

	const float DH = 2.1f;
	if (fabs(dis15 - 6.37f) < DH && fabs(dis14 - 5.18f) < DH &&
		fabs(dis25 - 5.18f) < DH && fabs(dis13 - 5.45f) < DH &&
		fabs(dis24 - 5.45f) < DH && fabs(dis35 - 5.45f) < DH)
		return 0;	// helix

	const float DS = 1.42f;
	if (fabs(dis15 - 13.0f) < DS && fabs(dis14 - 10.4f) < DS &&
		fabs(dis25 - 10.4f) < DS && fabs(dis13 - 6.1f) < DS &&
		fabs(dis24 - 6.1f) < DS && fabs(dis35 - 6.1f) < DS)
		return 1;	// strand

	if (dis15 < 8.2)
		return 2;	// turn

	return 3; // loop
	}

// 0=helix 1=strand 2=other
uint8_t chaq::get_ss3(const sid_t *distmx, uint M, uint L, uint pos)
	{
	uint8_t ss4 = get_ss4(distmx, M, L, pos);
	return ss4 <= 2 ? ss4 : 2;
	}

void chaq::get_ss4_str(const sid_t *distmx, uint M, uint L, string &ss)
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

void chaq::get_ss4_intseq(const sid_t *distmx, uint M, uint L, uint8_t *intseq)
	{
	for (uint pos = 0; pos < L; ++pos)
		intseq[pos] = get_ss4(distmx, M, L, pos);
	}

void chaq::create_nenvec(const sid_t* __restrict distmx, uint M, uint L,
	uint m, nnvec_t*& __restrict nnvec, sidvec_t*& __restrict nnsidvec)
	{
	nnvec = nnvec_t::newflat(L);
	nnsidvec = sidvec_t::newflat(L);
	uint16_t* __restrict v = nnvec->m_data;
	sid_t* __restrict dv = nnsidvec->m_data;

	for (uint i = 0; i < L; ++i)
		{
		v[i] = UINT16_MAX;
		dv[i] = UINT16_MAX;
		}

	for (uint i = 0; i < L; ++i)
		{
		uint dmax = L - 1 - i;
		if (dmax > M)
			dmax = M;

		for (uint d = m; d <= dmax; ++d)
			{
			uint j = i + d;
			sid_t sid = distmx[banded_ij_to_k(M, i, j)];

			if (sid < dv[i])
				{
				dv[i] = sid;
				v[i] = uint16_t(j);
				}

			if (sid < dv[j])
				{
				dv[j] = sid;
				v[j] = uint16_t(i);
				}
			}
		}
	}
