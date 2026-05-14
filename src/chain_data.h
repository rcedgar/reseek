#pragma once

#include "flat_dist_types.h"
#include "parasail.h"
#include "flat_chain.h"
#include "scratch_mem.h"

static const uint32_t bit_distmx =				(1 << 0);
static const uint32_t bit_mega_prof =			(1 << 1);
static const uint32_t bit_mega_prof_rev =		(1 << 2);
static const uint32_t bit_mega_pssm =			(1 << 3);
static const uint32_t bit_mega_pssm_rev =		(1 << 4);
static const uint32_t bit_parasail_prof =		(1 << 5);
static const uint32_t bit_parasail_prof_rev =	(1 << 6);
static const uint32_t bit_nu_codeseq =			(1 << 7);
static const uint32_t bit_nu_codeseq_rev =		(1 << 8);

static const uint32_t bits_query =
	bit_mega_pssm |
	bit_mega_prof |
	bit_parasail_prof |
	bit_parasail_prof_rev |
	bit_nu_codeseq |
	bit_nu_codeseq_rev;

static const uint32_t bits_target =
	bit_mega_prof |
	bit_nu_codeseq |
	bit_nu_codeseq_rev;

class chain_data
	{
public:
	static const uint32_t m_maxL;

public:
	string m_label;
	const flat_chain_t *m_chain = 0;
	uint m_L = 0;
	sid_t *m_distmx = 0;
	uint8_t *m_codeseq_nu = 0;
	uint8_t *m_codeseq_nu_rev = 0;
	uint8_t *m_mega_prof = 0;
	uint8_t *m_mega_prof_rev = 0;
	float *m_mega_pssm = 0;
	float *m_mega_pssm_rev = 0;
	parasail_profile_t *m_parasail_prof = 0;
	parasail_profile_t *m_parasail_prof_rev = 0;

public:
	static chain_data *from_chain(
		const flat_chain_t &chain,
		uint32_t bits,
		scratch_mem &mem,
		scratch_mem &scratch);

	static void fill_chain_data_vec(
		const vector<flat_chain_t *> &chains,
		uint32_t bits,
		chain_data **cdvec);

	static void get_object_counts(
		uint32_t bits,
		uint &n_distmx,
		uint &n_mega_prof,
		uint &n_mega_pssm,
		uint &n_codeseq);

	static void get_from_chain_bytes_per_pos(
		uint32_t bits,
		size_t &mem_bytes_per_pos,
		size_t &scratch_bytes_per_pos);

	//static size_t get_make_mega_prof_scratch_bytes_per_pos(uint32_t bits);

	static void make_mega_prof(
		const flat_chain_t &chain,
		const sid_t *distmx,
		uint8_t *mega_prof,
		size_t bytes,
		scratch_mem &mem,
		scratch_mem &scratch);

	static void log_mem_stats(chain_data **cdvec, uint n);
	};
