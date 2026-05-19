#include "myutils.h"
#include "parasail_nomalloc.h"

static size_t nomalloc_align_up(size_t n, size_t alignment)
	{
	asserta(alignment > 0);
	const size_t mask = alignment - 1;
	asserta((alignment & mask) == 0);
	return (n + mask) & ~mask;
	}

static __m256i *nomalloc_bump___m256i(char *&p, char *end, size_t alignment, int segLen)
	{
	p = (char *) nomalloc_align_up((size_t) p, alignment);
	const size_t nbytes = (size_t) segLen * sizeof(__m256i);
	asserta((size_t) (end - p) >= nbytes);
	__m256i *r = (__m256i *) p;
	p += nbytes;
	return r;
	}

size_t parasail_nomalloc_sw_striped_profile_avx2_256_16_workspace_bytes(int s1Len)
	{
	asserta(s1Len > 0);
	const int segWidth = 16;
	const int segLen = (s1Len + segWidth - 1) / segWidth;
	const size_t slab = (size_t) segLen * sizeof(__m256i);
	size_t total = 0;
	for (int i = 0; i < 4; ++i)
		{
		total = nomalloc_align_up(total, 32);
		total += slab;
		}
	return total;
	}

// Derived from parasail_sw_striped_profile_avx2_256_16 in parasail.cpp / sw_striped_avx2_256_16.c
int parasail_sw_striped_profile_avx2_256_16_nomalloc(
	const parasail_profile_t * const restrict profile,
	const char * const restrict s2, const int s2Len,
	const int open, const int gap,
	void *workspace, size_t workspace_bytes,
	int *out_end_query, int *out_end_ref,
	int *out_saturated)
	{
	int32_t i = 0;
	int32_t j = 0;
	int32_t k = 0;
	int32_t end_query = 0;
	int32_t end_ref = 0;
	int32_t s1Len = 0;
	const parasail_matrix_t *matrix = NULL;
	int32_t segWidth = 0;
	int32_t segLen = 0;
	__m256i *restrict vProfile = NULL;
	__m256i *restrict pvHStore = NULL;
	__m256i *restrict pvHLoad = NULL;
	__m256i *restrict pvHMax = NULL;
	__m256i *restrict pvE = NULL;
	__m256i vGapO;
	__m256i vGapE;
	__m256i vZero;
	int16_t bias = 0;
	int16_t score = 0;
	__m256i vBias;
	__m256i vMaxH;
	__m256i vMaxHUnit;
	int16_t maxp = 0;
	__m256i insert_mask;
	int saturated = 0;

	PARASAIL_CHECK_NULL(profile);
	PARASAIL_CHECK_NULL(profile->profile16.score);
	PARASAIL_CHECK_NULL(profile->matrix);
	PARASAIL_CHECK_GT0(profile->s1Len);
	PARASAIL_CHECK_NULL(s2);
	PARASAIL_CHECK_GT0(s2Len);
	PARASAIL_CHECK_GE0(open);
	PARASAIL_CHECK_GE0(gap);
	PARASAIL_CHECK_NULL(workspace);

	s1Len = profile->s1Len;
	segWidth = 16;
	segLen = (s1Len + segWidth - 1) / segWidth;

	const size_t need = parasail_nomalloc_sw_striped_profile_avx2_256_16_workspace_bytes(s1Len);
	asserta(workspace_bytes >= need);

	char *wp = (char *) workspace;
	char *wend = wp + workspace_bytes;
	pvHStore = nomalloc_bump___m256i(wp, wend, 32, segLen);
	pvHLoad = nomalloc_bump___m256i(wp, wend, 32, segLen);
	pvHMax = nomalloc_bump___m256i(wp, wend, 32, segLen);
	pvE = nomalloc_bump___m256i(wp, wend, 32, segLen);

	matrix = profile->matrix;
	vProfile = (__m256i *) profile->profile16.score;
	vGapO = _mm256_set1_epi16(open);
	vGapE = _mm256_set1_epi16(gap);
	vZero = _mm256_set1_epi16(0);
	bias = INT16_MIN;
	score = bias;
	vBias = _mm256_set1_epi16(bias);
	vMaxH = vBias;
	vMaxHUnit = vBias;
	maxp = INT16_MAX - (int16_t) (matrix->max + 1);
	insert_mask = _mm256_cmpgt_epi16(
		_mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
		vZero);

	parasail_memset___m256i(pvHStore, vBias, segLen);
	parasail_memset___m256i(pvE, vBias, segLen);

	for (j = 0; j < s2Len; ++j)
		{
		__m256i vE;
		__m256i vF;
		__m256i vH;
		const __m256i *vP = NULL;
		__m256i *pv = NULL;

		vF = vBias;

		vH = _mm256_slli_si256_rpl(pvHStore[segLen - 1], 2);
		vH = _mm256_blendv_epi8(vH, vBias, insert_mask);

		vP = vProfile + matrix->mapper[(unsigned char) s2[j]] * segLen;

		if (end_ref == j - 2)
			{
			pv = pvHMax;
			pvHMax = pvHLoad;
			pvHLoad = pvHStore;
			pvHStore = pv;
			}
		else
			{
			pv = pvHLoad;
			pvHLoad = pvHStore;
			pvHStore = pv;
			}

		for (i = 0; i < segLen; ++i)
			{
			vH = _mm256_adds_epi16(vH, _mm256_load_si256(vP + i));
			vE = _mm256_load_si256(pvE + i);

			vH = _mm256_max_epi16(vH, vE);
			vH = _mm256_max_epi16(vH, vF);
			_mm256_store_si256(pvHStore + i, vH);
			vMaxH = _mm256_max_epi16(vH, vMaxH);

			vH = _mm256_subs_epi16(vH, vGapO);
			vE = _mm256_subs_epi16(vE, vGapE);
			vE = _mm256_max_epi16(vE, vH);
			_mm256_store_si256(pvE + i, vE);

			vF = _mm256_subs_epi16(vF, vGapE);
			vF = _mm256_max_epi16(vF, vH);

			vH = _mm256_load_si256(pvHLoad + i);
			}

		for (k = 0; k < segWidth; ++k)
			{
			vF = _mm256_slli_si256_rpl(vF, 2);
			vF = _mm256_blendv_epi8(vF, vBias, insert_mask);
			for (i = 0; i < segLen; ++i)
				{
				vH = _mm256_load_si256(pvHStore + i);
				vH = _mm256_max_epi16(vH, vF);
				_mm256_store_si256(pvHStore + i, vH);
				vMaxH = _mm256_max_epi16(vH, vMaxH);
				vH = _mm256_subs_epi16(vH, vGapO);
				vF = _mm256_subs_epi16(vF, vGapE);
				if (!_mm256_movemask_epi8(_mm256_cmpgt_epi16(vF, vH)))
					goto end;
				}
			}
end:
		{
		__m256i vCompare = _mm256_cmpgt_epi16(vMaxH, vMaxHUnit);
		if (_mm256_movemask_epi8(vCompare))
			{
			score = _mm256_hmax_epi16_rpl(vMaxH);
			if (score > maxp)
				{
				saturated = 1;
				break;
				}
			vMaxHUnit = _mm256_set1_epi16(score);
			end_ref = j;
			}
		}
		}

	if (score == INT16_MAX)
		saturated = 1;

	if (saturated)
		{
		score = INT16_MAX;
		end_query = 0;
		end_ref = 0;
		}
	else
		{
		if (end_ref == j - 1)
			{
			__m256i *pv = pvHMax;
			pvHMax = pvHStore;
			pvHStore = pv;
			}
		else if (end_ref == j - 2)
			{
			__m256i *pv = pvHMax;
			pvHMax = pvHLoad;
			pvHLoad = pv;
			}
		{
		int16_t *t = (int16_t *) pvHMax;
		int32_t column_len = segLen * segWidth;
		end_query = s1Len - 1;
		for (i = 0; i < column_len; ++i, ++t)
			{
			if (*t == score)
				{
				int32_t temp = i / segWidth + i % segWidth * segLen;
				if (temp < end_query)
					end_query = temp;
				}
			}
		}
		}

	if (out_end_query != 0)
		*out_end_query = end_query;
	if (out_end_ref != 0)
		*out_end_ref = end_ref;
	if (out_saturated != 0)
		*out_saturated = saturated;

	return (int) score - (int) bias;
	}
