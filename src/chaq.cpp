#include "myutils.h"
#include "flat_base.h"
#include "chaq.h"
#include "flat_distmx.h"

static uint8_t get_aa4(char c)
	{
	c = toupper(c);
	if (c == 'G')
		return 0;
	if (strchr("AHPST", c) != 0)
		return 1;
	if (strchr("DEKNQR", c) != 0)
		return 2;
	return 3;
	}

static uint8_t get_aa3(char c)
	{
	c = toupper(c);
	if (c == 'G')
		return 0;
	if (strchr("CFILMVWY", c) != 0)
		return 1;
	return 2;
	}

void chaq::fill_distmx(
	cp_ic_t xyz,
	uint L,
	uint M,
	uint16_t *distmx)
	{
	fill_flat_distmx(xyz, L, M, distmx);
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
uint8_t chaq::get_ss3(const sid_t * __restrict distmx, uint M, uint L, uint pos)
	{
	uint8_t ss4 = get_ss4(distmx, M, L, pos);
	return ss4 <= 2 ? ss4 : 2;
	}

void chaq::get_ss4_str(const sid_t * __restrict distmx, uint M, uint L, string &ss)
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

void chaq::get_ss4_codeseq(cp_sid_t distmx, uint M, uint L, p_uint8_t intseq)
	{
	for (uint pos = 0; pos < L; ++pos)
		intseq[pos] = get_ss4(distmx, M, L, pos);
	}

void chaq::fill_nenvec(
	cp_sid_t distmx,
	uint L,
	uint M,
	uint m,
	p_uint16_t nen,
	p_uint16_t nensid)
	{
	for (uint i = 0; i < L; ++i)
		{
		nen[i] = UINT16_MAX;
		nensid[i] = UINT16_MAX;
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

			if (sid < nensid[i])
				{
				nensid[i] = sid;
				nen[i] = uint16_t(j);
				}

			if (sid < nensid[j])
				{
				nensid[j] = sid;
				nen[j] = uint16_t(i);
				}
			}
		}
	}

void chaq::fill_nen_ren_vecs(
	cp_uint16_t pens,
	cp_uint16_t mens,
	cp_sid_t pensids,
	cp_sid_t mensids,
	uint L,
	p_uint16_t nens,
	p_uint16_t rens,
	p_uint16_t nensids,
	p_uint16_t rensids)
	{
	for (uint i = 0; i < L; ++i)
		{
		sid_t pensid = pensids[i];
		sid_t mensid = mensids[i];
		if (pensid <= mensid)
			{
			nens[i] = pens[i];
			rens[i] = mens[i];
			nensids[i] = pensids[i];
			rensids[i] = mensids[i];
			}
		else
			{
			nens[i] = mens[i];
			rens[i] = pens[i];
			nensids[i] = mensids[i];
			rensids[i] = pensids[i];
			}
		}
	}

void chaq::fill_pen_men_vecs(
	cp_sid_t distmx,
	uint L,
	uint M,
	uint m,
	p_uint16_t pen,
	p_uint16_t pensid,
	p_uint16_t men,
	p_uint16_t mensid)
	{
	for (uint i = 0; i < L; ++i)
		{
		pen[i]  = UINT16_MAX;
		pensid[i] = UINT16_MAX;
		men[i]  = UINT16_MAX;
		mensid[i] = UINT16_MAX;
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

			// forward for i: j > i
			if (sid < pensid[i])
				{
				pensid[i] = sid;
				pen[i] = uint16_t(j);
				}

			// reverse for j: i < j
			if (sid < mensid[j])
				{
				mensid[j] = sid;
				men[j] = uint16_t(i);
				}
			}
		}
	}

void chaq::get_pm_codeseq(cp_sid_t pensids, cp_sid_t mensids, uint L, p_uint8_t codeseq)
	{
	for (uint i = 0; i < L; ++i)
		codeseq[i] = (pensids[i] <= mensids[i] ? 0 : 1);
	}
