#pragma once

#include "flat_base.h"
#include "flat_chain.h"
#include "fan.h"

/***
Chain quantizer / quantifier
All member functions are static.
Function arguments are pointers to data, not to flat_base objects.
Consuming code responsible for memory ownership and lifefime.
***/
class chaq
	{
private:
	chaq() = delete;

public:
	static void fill_distmx(
		cp_ic_t xyz,
		uint L,
		uint M,
		uint16_t *distmx);

	static void fill_nenvec(
		cp_sid_t distmx,
		uint L,
		uint M,
		uint m,
		p_uint16_t nnvec,
		p_uint16_t nnsidvec);

	static void fill_nen_ren_vecs(
		cp_uint16_t pens,
		cp_uint16_t mens,
		cp_sid_t pensids,
		cp_sid_t mensids,
		uint L,
		p_uint16_t nens,
		p_uint16_t rens,
		p_uint16_t nensids,
		p_uint16_t rensids);

	static void fill_pen_men_vecs(
		cp_sid_t distmx,
		uint L,
		uint M,
		uint m,
		p_uint16_t penvec,
		p_uint16_t pensidvec,
		p_uint16_t menvec,
		p_uint16_t mensidvec);

	static uint8_t get_ss3(cp_sid_t distmx, uint M, uint L, uint pos);
	static uint8_t get_ss4(cp_sid_t distmx, uint M, uint L, uint pos);
	static void get_ss4_str(cp_sid_t distmx, uint M, uint L, string &ss);

	static void get_ss3_codeseq(cp_sid_t distmx, uint M, uint L, p_uint8_t codeseq);
	static void get_ss4_codeseq(cp_sid_t distmx, uint M, uint L, p_uint8_t codeseq);
	};
