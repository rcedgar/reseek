#pragma once

#include "flat_base.h"
#include "flat_chain.h"

/***
Chain quantizer / quantifier
All member functions are static.
Function arguments are myptr's.
Consuming code responsible for ownership and lifefime.
***/
class chaq
	{
private:
	chaq() = delete;

public:
	static void create_distmx(const flat_chain_t *chain, 
		chaindistmx_t*& dm, uint M);
	static void create_nenvec(const sid_t *distmx, uint M, uint L,
		uint m, nnvec_t*& nnvec, sidvec_t*& nnsidvec);
	static uint8_t get_ss3(const sid_t *distmx, uint M, uint L, uint pos);
	static uint8_t get_ss4(const sid_t *distmx, uint M, uint L, uint pos);
	static void get_ss4_str(const sid_t *distmx, uint M, uint L, string &ss);
	static void get_ss4_intseq(const sid_t *distmx, uint M, uint L, uint8_t *intseq);
	};
