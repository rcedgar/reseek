#pragma once

#include "parasail_nomalloc.h"

struct query_data
	{
	const flat_chain_t *m_chain = 0;
	sid_t *m_distmx = 0;
	uint8_t *m_codeseq_nu = 0;
	uint8_t *m_codeseq_kappa = 0;
	uint8_t *m_codeseq_nu_rev = 0;
	uint8_t *m_mega_prof = 0;
	float *m_mega_pssm = 0;
	float *m_mega_pssm_rev = 0;
	parasail_profile_t *m_parasail_prof = 0;
	parasail_profile_t *m_parasail_prof_rev = 0;
	};

struct db_data
	{
	const flat_chain_t *m_chain = 0;
	sid_t *m_distmx = 0;
	uint8_t *m_codeseq_nu = 0;
	uint8_t *m_codeseq_nu_rev = 0;
	parasail_profile_t *m_parasail_prof = 0;

	static void delete_db_data(db_data *dd);
	};
