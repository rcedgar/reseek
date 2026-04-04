#pragma once

#include "flat_base.h"
#include "flat_chain.h"
#include "fan.h"
#include "alpha.h"

class sec_kmeans;

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
	static bool feature_is_binned(FAN fan)
		{
		switch (fan)
			{
		case FAN_nendist:
		case FAN_rendist:
		case FAN_pendist:
		case FAN_mendist:
		case FAN_pack:
		case FAN_ppack:
		case FAN_mpack:
			return true;
			}
		return false;
		}

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

	static void get_aa3_codeseq(const char *aacharseq, uint L, p_uint8_t codeseq);
	static void get_aa4_codeseq(const char *aacharseq, uint L, p_uint8_t codeseq);

	static void get_pm_codeseq(cp_sid_t pensids, cp_sid_t mensids, uint L, p_uint8_t codeseq);

	static void get_packing_values(cp_sid_t distmx, uint M, uint L,
		uint maxsid, bool include_plus, bool include_minus,
		p_uint16_t values);

	static void get_packing_codeseq(cp_sid_t distmx, uint M, uint L, 
		uint maxsid, bool include_plus, bool include_minus, p_uint8_t codeseq);

	static void slow_get_values(
		const flat_chain_t *chain,
		FAN fan,
		uint alpha_size,
		uint M,
		uint m,
		p_uint16_t values);

	static void slow_get_codeseq_binned(
		const flat_chain_t *chain,
		FAN fan,
		uint alpha_size,
		uint M,
		uint m,
		p_uint8_t codeseq);

	static void slow_get_codeseq_discrete(
		const flat_chain_t *chain,
		FAN fan,
		uint alpha_size,
		uint M,
		uint m,
		uint8_t undef_code,
		p_uint8_t codeseq);

	static void slow_get_charseq_binned(
		const flat_chain_t *chain,
		FAN fan,
		uint8_t alpha_size,
		uint M,
		uint m,
		cp_uint16_t thresholds,
		uint16_t undef_value,
		char *charseq);

	static void slow_get_charseq_discrete(
		const flat_chain_t *chain,
		FAN fan,
		uint8_t alpha_size,
		uint M,
		uint m,
		uint8_t undef_code,
		char *charseq);

	static void get_sec_codeseq(
		uint alpha_size,
		cp_sid_t distmx,
		uint M,
		uint L,
		p_uint8_t codeseq);

	static const uint8_t *get_char2letter(uint alpha_size)
		{
		return (alpha_size == 20 ? g_CharToLetterAmino : g_CharToLetterMu);
		}

	static const unsigned char *get_letter2char(uint alpha_size)
		{
		return (alpha_size == 20 ? g_LetterToCharAmino : g_LetterToCharMu);
		}

	static void codeseq2charseq(
		const uint8_t *codeseq,
		uint L,
		uint alpha_size,
		char *charseq);

	static void charseq2codeseq(
		const char *charseq,
		uint L,
		uint alpha_size,
		uint8_t *codeseq);

	static cp_uint16_t get_thresholds(FAN fan, uint alpha_size);
	static uint16_t get_undef_value(FAN fan, uint alpha_size);
	static uint8_t get_undef_code(FAN fan, uint alpha_size);
	};
