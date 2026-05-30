#include "myutils.h"
#include "nu_filter.h"
#include "flat_nu_aligner.h"

uint nu_filter::m_nchain = 0;
const vector<string> *nu_filter::m_ptr_labels = 0;
parasail_profile_t **nu_filter::m_parasail_profs = 0;
parasail_profile_t **nu_filter::m_parasail_prof_revs = 0;
const uint *nu_filter::m_lengths = 0;

void nu_filter::set_query_parasail_profiles(
	const vector<string> &labels,
	const uint8_t **codeseqs_nu,
	const uint *lengths,
	uint nchain)
	{
	asserta(m_nchain == 0);
	asserta(labels.size() == nchain);
	m_nchain = nchain;
	m_ptr_labels = &labels;
	m_lengths = lengths;
	m_parasail_profs = myalloc(parasail_profile_t *, nchain);
	m_parasail_prof_revs = myalloc(parasail_profile_t *, nchain);

	for (uint chainidx = 0; chainidx < nchain; ++chainidx)
		{
		const uint8_t *nu_codeseq = codeseqs_nu[chainidx];
		uint L = lengths[chainidx];
		parasail_profile_t *prof = parasail_profile_create_avx_256_16(
			(const char *) nu_codeseq, L, &flat_nu_aligner::m_matrix);
		uint8_t *revseq = myalloc(uint8_t, L);
		for (uint i = 0; i < L; ++i)
			revseq[i] = nu_codeseq[L-i-1];
		parasail_profile_t *prof_rev = parasail_profile_create_avx_256_16(
			(const char *) nu_codeseq, L, &flat_nu_aligner::m_matrix);
		myfree(revseq);
		} 
	}
