#pragma once

#include "flat_dist_types.h"
#include "parasail.h"
#include "flat_chain.h"
#include "tdata.h"

class qdata : public tdata
	{
public:
	static qdata *m_qdatavec;
	static uint32_t m_nq;

public:
	uint8_t *m_codeseq_nu_rev = 0;
	float *m_mega_pssm = 0;
	parasail_profile_t *m_parasail_prof = 0;
	parasail_profile_t *m_parasail_prof_rev = 0;

public:
	qdata()
		{
		m_isq = true;
		}
	};
