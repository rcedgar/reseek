#pragma once

#include "parasail_nomalloc.h"

class nu_filter
	{
public:
	static uint m_nchain;
	static const vector<string> *m_ptr_labels;
	static parasail_profile_t **m_parasail_profs;
	static parasail_profile_t **m_parasail_prof_revs;
	static const uint *m_lengths;

public:
	static void set_query_parasail_profiles(
		const vector<string> &labels,
		const uint8_t **codeseqs_nu,
		const uint *lengths,
		uint nchain);
	};