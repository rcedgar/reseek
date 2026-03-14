#pragma once

#include "flat_base.h"
#include "flat_chain.h"

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
		const ic_t *xyz,
		uint L,
		uint M,
		uint16_t *distmx);

	static void fill_nenvec(
		const sid_t* __restrict distmx,
		uint L,
		uint M,
		uint m,
		uint16_t* __restrict nnvec,
		uint16_t* __restrict nnsidvec);

	static void fill_nen_pen_vecs(
		const sid_t* __restrict distmx,
		uint L,
		uint M,
		uint m,
		uint16_t* __restrict nnvec,
		uint16_t* __restrict nnsidvec,
		uint16_t* __restrict renvec,
		uint16_t* __restrict renidvec);

	static uint8_t get_ss3(const sid_t *distmx, uint M, uint L, uint pos);
	static uint8_t get_ss4(const sid_t *distmx, uint M, uint L, uint pos);
	static void get_ss4_str(const sid_t *distmx, uint M, uint L, string &ss);
	static void get_ss4_intseq(const sid_t *distmx, uint M, uint L, uint8_t *intseq);
	};
