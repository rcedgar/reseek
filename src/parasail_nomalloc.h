#pragma once

#include "parasail.h"
#include <stddef.h>

/* Bytes required for workspace of parasail_sw_striped_profile_avx2_256_16_nomalloc.
   s1Len is query length (profile->s1Len). Internal alignment padding included. */
uint parasail_nomalloc_sw_striped_profile_avx2_256_16_workspace_bytes(int s1Len);

/* Score-only striped SW (16-bit) using caller-provided workspace (no malloc/free).
   Returns alignment score (same as parasail_result_t::score from the malloc version).
   If out_saturated is non-null, set to 1 when the malloc version would set PARASAIL_FLAG_SATURATED.
   out_end_query / out_end_ref are optional. */
int parasail_sw_striped_profile_avx2_256_16_nomalloc(
	const parasail_profile_t * const restrict profile,
	const char * const restrict s2, const int s2Len,
	const int open, const int gap,
	void *workspace, size_t workspace_bytes,
	int *out_end_query = 0, int *out_end_ref = 0,
	int *out_saturated = 0);
