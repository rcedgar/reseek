#include "myutils.h"
#include "parasail.h"
#include "timing.h"

/***
This source code is partly derived from the parasail library
https://github.com/jeffdaily/parasail
Author jeffrey.daily@gmail.com
Copyright (c) 2015 Battelle Memorial Institute.
***/

parasail_profile_t* parasail_profile_new(
        const char * s1, const int s1Len, const parasail_matrix_t *matrix)
{
    /* declare all variables */
    parasail_profile_t *profile = NULL;

    asserta(matrix != 0);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        asserta(s1 != 0);
    }

    PARASAIL_NEW(profile, parasail_profile_t);

    profile->s1 = s1;
    profile->s1Len = s1Len;
    profile->matrix = matrix;
    profile->profile8.score = NULL;
    profile->profile8.matches = NULL;
    profile->profile8.similar = NULL;
    profile->profile16.score = NULL;
    profile->profile16.matches = NULL;
    profile->profile16.similar = NULL;
    profile->profile32.score = NULL;
    profile->profile32.matches = NULL;
    profile->profile32.similar = NULL;
    profile->profile64.score = NULL;
    profile->profile64.matches = NULL;
    profile->profile64.similar = NULL;
    profile->free = NULL;
    profile->stop = INT32_MAX;

    if (matrix->type == PARASAIL_MATRIX_TYPE_PSSM) {
        profile->s1Len = matrix->length;
    }

    return profile;
}

void* parasail_memalign(size_t alignment, size_t size)
{
    void *ptr = NULL;
#ifdef _MSC_VER
    ptr = _aligned_malloc(size, alignment);
#else
    ptr = aligned_alloc(alignment, size);
#endif
	asserta(ptr != 0);
    return ptr;
}

void parasail_free(void *ptr)
{
#ifdef _MSC_VER
	_aligned_free(ptr);
#else
	free(ptr);
#endif
}

parasail_profile_t * parasail_profile_create_avx_256_8(
        const char * const restrict s1, const int _s1Len,
        const parasail_matrix_t *matrix)
{
    int s1Len = 0;
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    int32_t n = 0; /* number of amino acids in table */
    const int32_t segWidth = 32; /* number of values in vector unit */
    int32_t segLen = 0;
    __m256i* restrict vProfile = NULL;
    int32_t index = 0;
    parasail_profile_t *profile = NULL;

    asserta(matrix != 0);
    /* s1 is ignored for pssm, required for square */
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        asserta(s1 != 0);
    }

    s1Len = matrix->type == PARASAIL_MATRIX_TYPE_SQUARE ? _s1Len : matrix->length;
    n = matrix->size;
    segLen = (s1Len + segWidth - 1) / segWidth;
    vProfile = parasail_memalign___m256i(32, n * segLen);
    if (!vProfile) return NULL;
    profile = parasail_profile_new(s1, s1Len, matrix);
    if (!profile) return NULL;

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            __m256i_8_t t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
                    t.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                }
                else {
                    t.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*j+matrix->mapper[(unsigned char)matrix->alphabet[k]]];
                }
                j += segLen;
            }
            _mm256_store_si256(&vProfile[index], t.m);
            ++index;
        }
    }

    profile->profile8.score = vProfile;
    profile->free = &parasail_free___m256i;
    return profile;
}

parasail_result_t* parasail_result_new()
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    PARASAIL_NEW(result, parasail_result_t);

    result->score = 0;
    result->end_query = 0;
    result->end_ref = 0;
    result->flag = 0;

    return result;
}

parasail_result_t* parasail_result_new_trace(const int a, const int b, const size_t alignment, const size_t size)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    asserta(a > 0 && b > 0);

    /* allocate struct to hold memory */
    result = parasail_result_new();
    if (!result) return NULL;

    PARASAIL_NEW(result->trace, parasail_result_extra_trace_t);
    result->trace->trace_table = parasail_memalign(alignment, size*a*b);
    if (!result->trace->trace_table) return NULL;
    result->trace->trace_ins_table = NULL;
    result->trace->trace_del_table = NULL;

    return result;
}

parasail_result_t* parasail_sw_trace_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    /* declare local variables */
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int32_t s1Len = 0;
    const parasail_matrix_t *matrix = NULL;
    int32_t segWidth = 0;
    int32_t segLen = 0;
    __m256i* restrict vProfile = NULL;
    __m256i* restrict pvHStore = NULL;
    __m256i* restrict pvHLoad = NULL;
    __m256i* restrict pvE = NULL;
    __m256i* restrict pvEaStore = NULL;
    __m256i* restrict pvEaLoad = NULL;
    __m256i* restrict pvHT = NULL;
    __m256i* restrict pvHMax = NULL;
    __m256i vGapO;
    __m256i vGapE;
    __m256i vZero;
    int8_t score = 0;
    __m256i vMaxH;
    __m256i vMaxHUnit;
    int8_t maxp = 0;
    parasail_result_t *result = NULL;
    __m256i vTZero;
    __m256i vTIns;
    __m256i vTDel;
    __m256i vTDiag;
    __m256i vTDiagE;
    __m256i vTInsE;
    __m256i vTDiagF;
    __m256i vTDelF;
    __m256i vTMask;
    __m256i vFTMask;

    /* validate inputs */
    PARASAIL_CHECK_NULL(profile);
    PARASAIL_CHECK_NULL(profile->profile8.score);
    PARASAIL_CHECK_NULL(profile->matrix);
    PARASAIL_CHECK_GT0(profile->s1Len);
    PARASAIL_CHECK_NULL(s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);

    /* initialize stack variables */
    i = 0;
    j = 0;
    k = 0;
    end_query = 0;
    end_ref = 0;
    s1Len = profile->s1Len;
    matrix = profile->matrix;
    segWidth = 32; /* number of values in vector unit */
    segLen = (s1Len + segWidth - 1) / segWidth;
    vProfile = (__m256i*)profile->profile8.score;
    vGapO = _mm256_set1_epi8(open);
    vGapE = _mm256_set1_epi8(gap);
    vZero = _mm256_setzero_si256();
    score = INT8_MIN;
    vMaxH = vZero;
    vMaxHUnit = vZero;
    maxp = INT8_MAX - (int8_t)(matrix->max+1);
    vTZero = _mm256_set1_epi8(PARASAIL_ZERO);
    vTIns  = _mm256_set1_epi8(PARASAIL_INS);
    vTDel  = _mm256_set1_epi8(PARASAIL_DEL);
    vTDiag = _mm256_set1_epi8(PARASAIL_DIAG);
    vTDiagE= _mm256_set1_epi8(PARASAIL_DIAG_E);
    vTInsE = _mm256_set1_epi8(PARASAIL_INS_E);
    vTDiagF= _mm256_set1_epi8(PARASAIL_DIAG_F);
    vTDelF = _mm256_set1_epi8(PARASAIL_DEL_F);
    vTMask = _mm256_set1_epi8(PARASAIL_ZERO_MASK);
    vFTMask= _mm256_set1_epi8(PARASAIL_F_MASK);

    /* initialize result */
    result = parasail_result_new_trace(segLen, s2Len, 32, sizeof(__m256i));
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_32;

    /* initialize heap variables */
    pvHStore = parasail_memalign___m256i(32, segLen);
    pvHLoad =  parasail_memalign___m256i(32, segLen);
    pvE = parasail_memalign___m256i(32, segLen);
    pvEaStore = parasail_memalign___m256i(32, segLen);
    pvEaLoad = parasail_memalign___m256i(32, segLen);
    pvHT = parasail_memalign___m256i(32, segLen);
    pvHMax = parasail_memalign___m256i(32, segLen);

    /* validate heap variables */
    if (!pvHStore) return NULL;
    if (!pvHLoad) return NULL;
    if (!pvE) return NULL;
    if (!pvEaStore) return NULL;
    if (!pvEaLoad) return NULL;
    if (!pvHT) return NULL;
    if (!pvHMax) return NULL;

    /* initialize H and E */
    parasail_memset___m256i(pvHStore, vZero, segLen);
    parasail_memset___m256i(pvE, _mm256_set1_epi8(-open), segLen);
    parasail_memset___m256i(pvEaStore, _mm256_set1_epi8(-open), segLen);

    for (i=0; i<segLen; ++i) {
        arr_store(result->trace->trace_table, vTDiagE, i, segLen, 0);
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vEF_opn;
        __m256i vE;
        __m256i vE_ext;
        __m256i vF;
        __m256i vF_ext;
        __m256i vFa;
        __m256i vFa_ext;
        __m256i vH;
        __m256i vH_dag;
        const __m256i* vP = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop. */
        vF = _mm256_subs_epi8(vZero,vGapO);

        /* load final segment of pvHStore and shift left by 1 bytes */
        vH = _mm256_load_si256(&pvHStore[segLen - 1]);
        vH = _mm256_slli_si256_rpl(vH, 1);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        if (end_ref == j-2) {
            /* Swap in the max buffer. */
            SWAP3(pvHMax,  pvHLoad,  pvHStore)
            SWAP(pvEaLoad,  pvEaStore)
        }
        else {
            /* Swap the 2 H buffers. */
            SWAP(pvHLoad,  pvHStore)
            SWAP(pvEaLoad,  pvEaStore)
        }

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vE = _mm256_load_si256(pvE + i);

            /* Get max from vH, vE and vF. */
            vH_dag = _mm256_adds_epi8(vH, _mm256_load_si256(vP + i));
            vH_dag = _mm256_max_epi8(vH_dag, vZero);
            vH = _mm256_max_epi8(vH_dag, vE);
            vH = _mm256_max_epi8(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);

            {
                __m256i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                __m256i cond_zero = _mm256_cmpeq_epi8(vH, vZero);
                __m256i case1 = _mm256_cmpeq_epi8(vH, vH_dag);
                __m256i case2 = _mm256_cmpeq_epi8(vH, vF);
                __m256i vT = _mm256_blendv_epi8(
                        _mm256_blendv_epi8(vTIns, vTDel, case2),
                        _mm256_blendv_epi8(vTDiag, vTZero, cond_zero),
                        case1);
                _mm256_store_si256(pvHT + i, vT);
                vT = _mm256_or_si256(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i, segLen, j);
            }
            vMaxH = _mm256_max_epi8(vH, vMaxH);
            vEF_opn = _mm256_subs_epi8(vH, vGapO);

            /* Update vE value. */
            vE_ext = _mm256_subs_epi8(vE, vGapE);
            vE = _mm256_max_epi8(vEF_opn, vE_ext);
            _mm256_store_si256(pvE + i, vE);
            {
                __m256i vEa = _mm256_load_si256(pvEaLoad + i);
                __m256i vEa_ext = _mm256_subs_epi8(vEa, vGapE);
                vEa = _mm256_max_epi8(vEF_opn, vEa_ext);
                _mm256_store_si256(pvEaStore + i, vEa);
                if (j+1<s2Len) {
                    __m256i cond = _mm256_cmpgt_epi8(vEF_opn, vEa_ext);
                    __m256i vT = _mm256_blendv_epi8(vTInsE, vTDiagE, cond);
                    arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                }
            }

            /* Update vF value. */
            vF_ext = _mm256_subs_epi8(vF, vGapE);
            vF = _mm256_max_epi8(vEF_opn, vF_ext);
            if (i+1<segLen) {
                __m256i vTAll = arr_load(result->trace->trace_table, i+1, segLen, j);
                __m256i cond = _mm256_cmpgt_epi8(vEF_opn, vF_ext);
                __m256i vT = _mm256_blendv_epi8(vTDelF, vTDiagF, cond);
                vT = _mm256_or_si256(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i+1, segLen, j);
            }

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        vFa_ext = vF_ext;
        vFa = vF;
        for (k=0; k<segWidth; ++k) {
            __m256i vHp = _mm256_load_si256(&pvHLoad[segLen - 1]);
            vHp = _mm256_slli_si256_rpl(vHp, 1);
            vEF_opn = _mm256_slli_si256_rpl(vEF_opn, 1);
            vEF_opn = _mm256_insert_epi8_rpl(vEF_opn, -open, 0);
            vF_ext = _mm256_slli_si256_rpl(vF_ext, 1);
            vF_ext = _mm256_insert_epi8_rpl(vF_ext, NEG_INF, 0);
            vF = _mm256_slli_si256_rpl(vF, 1);
            vF = _mm256_insert_epi8_rpl(vF, -open, 0);
            vFa_ext = _mm256_slli_si256_rpl(vFa_ext, 1);
            vFa_ext = _mm256_insert_epi8_rpl(vFa_ext, NEG_INF, 0);
            vFa = _mm256_slli_si256_rpl(vFa, 1);
            vFa = _mm256_insert_epi8_rpl(vFa, -open, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi8(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
                {
                    __m256i vTAll;
                    __m256i vT;
                    __m256i case1;
                    __m256i case2;
                    __m256i cond;
                    vHp = _mm256_adds_epi8(vHp, _mm256_load_si256(vP + i));
                    vHp = _mm256_max_epi8(vHp, vZero);
                    case1 = _mm256_cmpeq_epi8(vH, vHp);
                    case2 = _mm256_cmpeq_epi8(vH, vF);
                    cond = _mm256_andnot_si256(case1,case2);
                    vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    vT = _mm256_load_si256(pvHT + i);
                    vT = _mm256_blendv_epi8(vT, vTDel, cond);
                    _mm256_store_si256(pvHT + i, vT);
                    vTAll = _mm256_and_si256(vTAll, vTMask);
                    vTAll = _mm256_or_si256(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                vMaxH = _mm256_max_epi8(vH, vMaxH);
                /* Update vF value. */
                {
                    __m256i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    __m256i cond = _mm256_cmpgt_epi8(vEF_opn, vFa_ext);
                    __m256i vT = _mm256_blendv_epi8(vTDelF, vTDiagF, cond);
                    vTAll = _mm256_and_si256(vTAll, vFTMask);
                    vTAll = _mm256_or_si256(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                vEF_opn = _mm256_subs_epi8(vH, vGapO);
                vF_ext = _mm256_subs_epi8(vF, vGapE);
                {
                    __m256i vEa = _mm256_load_si256(pvEaLoad + i);
                    __m256i vEa_ext = _mm256_subs_epi8(vEa, vGapE);
                    vEa = _mm256_max_epi8(vEF_opn, vEa_ext);
                    _mm256_store_si256(pvEaStore + i, vEa);
                    if (j+1<s2Len) {
                        __m256i cond = _mm256_cmpgt_epi8(vEF_opn, vEa_ext);
                        __m256i vT = _mm256_blendv_epi8(vTInsE, vTDiagE, cond);
                        arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                    }
                }
                if (! _mm256_movemask_epi8(
                            _mm256_or_si256(
                                _mm256_cmpgt_epi8(vF_ext, vEF_opn),
                                _mm256_cmpeq_epi8(vF_ext, vEF_opn))))
                    goto end;
                /*vF = _mm256_max_epi8(vEF_opn, vF_ext);*/
                vF = vF_ext;
                vFa_ext = _mm256_subs_epi8(vFa, vGapE);
                vFa = _mm256_max_epi8(vEF_opn, vFa_ext);
                vHp = _mm256_load_si256(pvHLoad + i);
            }
        }
end:
        {
        }

        {
            __m256i vCompare = _mm256_cmpgt_epi8(vMaxH, vMaxHUnit);
            if (_mm256_movemask_epi8(vCompare)) {
                score = _mm256_hmax_epi8_rpl(vMaxH);
                /* if score has potential to overflow, abort early */
                if (score > maxp) {
                    result->flag |= PARASAIL_FLAG_SATURATED;
                    break;
                }
                vMaxHUnit = _mm256_set1_epi8(score);
                end_ref = j;
            }
        }

        /*if (score == stop) break;*/
    }

    if (score == INT8_MAX) {
        result->flag |= PARASAIL_FLAG_SATURATED;
    }

    if (parasail_result_is_saturated(result)) {
        score = 0;
        end_query = 0;
        end_ref = 0;
    }
    else {
        if (end_ref == j-1) {
            /* end_ref was the last store column */
            SWAP(pvHMax,  pvHStore)
        }
        else if (end_ref == j-2) {
            /* end_ref was the last load column */
            SWAP(pvHMax,  pvHLoad)
        }
        /* Trace the alignment ending position on read. */
        {
            int8_t *t = (int8_t*)pvHMax;
            int32_t column_len = segLen * segWidth;
            end_query = s1Len - 1;
            for (i = 0; i<column_len; ++i, ++t) {
                if (*t == score) {
                    int32_t temp = i / segWidth + i % segWidth * segLen;
                    if (temp < end_query) {
                        end_query = temp;
                    }
                }
            }
        }
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(pvHMax);
    parasail_free(pvHT);
    parasail_free(pvEaLoad);
    parasail_free(pvEaStore);
    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

parasail_result_t* parasail_sw_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    /* declare local variables */
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int32_t s1Len = 0;
    const parasail_matrix_t *matrix = NULL;
    int32_t segWidth = 0;
    int32_t segLen = 0;
#ifdef PARASAIL_ROWCOL
    int32_t offset = 0;
    int32_t position = 0;
#endif
    __m256i* restrict vProfile = NULL;
    __m256i* restrict pvHStore = NULL;
    __m256i* restrict pvHLoad = NULL;
    __m256i* restrict pvHMax = NULL;
    __m256i* restrict pvE = NULL;
    __m256i vGapO;
    __m256i vGapE;
    __m256i vZero;
    int8_t bias = 0;
    int8_t score = 0;
    __m256i vBias;
    __m256i vMaxH;
    __m256i vMaxHUnit;
    int8_t maxp = 0;
    __m256i insert_mask;
    /*int8_t stop = 0;*/
    parasail_result_t *result = NULL;

    /* validate inputs */
    PARASAIL_CHECK_NULL(profile);
    PARASAIL_CHECK_NULL(profile->profile8.score);
    PARASAIL_CHECK_NULL(profile->matrix);
    PARASAIL_CHECK_GT0(profile->s1Len);
    PARASAIL_CHECK_NULL(s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);

    /* initialize stack variables */
    i = 0;
    j = 0;
    k = 0;
    end_query = 0;
    end_ref = 0;
    s1Len = profile->s1Len;
    matrix = profile->matrix;
    segWidth = 32; /* number of values in vector unit */
    segLen = (s1Len + segWidth - 1) / segWidth;
#ifdef PARASAIL_ROWCOL
    offset = (s1Len - 1) % segLen;
    position = (segWidth - 1) - (s1Len - 1) / segLen;
#endif
    vProfile = (__m256i*)profile->profile8.score;
    vGapO = _mm256_set1_epi8(open);
    vGapE = _mm256_set1_epi8(gap);
    vZero = _mm256_set1_epi8(0);
    bias = INT8_MIN;
    score = bias;
    vBias = _mm256_set1_epi8(bias);
    vMaxH = vBias;
    vMaxHUnit = vBias;
    maxp = INT8_MAX - (int8_t)(matrix->max+1);
    insert_mask = _mm256_cmpgt_epi8(
            _mm256_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
            vZero);
    /*stop = profile->stop == INT32_MAX ?  INT8_MAX : (int8_t)profile->stop-bias;*/

    /* initialize result */
#ifdef PARASAIL_TABLE
    result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    result = parasail_result_new_rowcol1(segLen*segWidth, s2Len);
#else
    result = parasail_result_new();
#endif
#endif
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_32;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    pvHStore = parasail_memalign___m256i(32, segLen);
    pvHLoad = parasail_memalign___m256i(32, segLen);
    pvHMax = parasail_memalign___m256i(32, segLen);
    pvE = parasail_memalign___m256i(32, segLen);

    /* validate heap variables */
    if (!pvHStore) return NULL;
    if (!pvHLoad) return NULL;
    if (!pvHMax) return NULL;
    if (!pvE) return NULL;

    /* initialize H and E */
    parasail_memset___m256i(pvHStore, vBias, segLen);
    parasail_memset___m256i(pvE, vBias, segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vE;
        __m256i vF;
        __m256i vH;
        const __m256i* vP = NULL;
        __m256i* pv = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        vF = vBias;

        /* load final segment of pvHStore and shift left by 1 bytes */
        vH = _mm256_slli_si256_rpl(pvHStore[segLen - 1], 1);
        vH = _mm256_blendv_epi8(vH, vBias, insert_mask);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        if (end_ref == j-2) {
            /* Swap in the max buffer. */
            pv = pvHMax;
            pvHMax = pvHLoad;
            pvHLoad = pvHStore;
            pvHStore = pv;
        }
        else {
            /* Swap the 2 H buffers. */
            pv = pvHLoad;
            pvHLoad = pvHStore;
            pvHStore = pv;
        }

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = _mm256_adds_epi8(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi8(vH, vE);
            vH = _mm256_max_epi8(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->tables->score_table, vH, i, segLen, j, s2Len, bias);
#endif
            vMaxH = _mm256_max_epi8(vH, vMaxH);

            /* Update vE value. */
            vH = _mm256_subs_epi8(vH, vGapO);
            vE = _mm256_subs_epi8(vE, vGapE);
            vE = _mm256_max_epi8(vE, vH);
            _mm256_store_si256(pvE + i, vE);

            /* Update vF value. */
            vF = _mm256_subs_epi8(vF, vGapE);
            vF = _mm256_max_epi8(vF, vH);

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            vF = _mm256_slli_si256_rpl(vF, 1);
            vF = _mm256_blendv_epi8(vF, vBias, insert_mask);
            for (i=0; i<segLen; ++i) {
                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi8(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si256(result->tables->score_table, vH, i, segLen, j, s2Len, bias);
#endif
                vMaxH = _mm256_max_epi8(vH, vMaxH);
                vH = _mm256_subs_epi8(vH, vGapO);
                vF = _mm256_subs_epi8(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi8(vF, vH))) goto end;
                /*vF = _mm256_max_epi8(vF, vH);*/
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm256_load_si256(pvHStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl(vH, 1);
            }
            result->rowcols->score_row[j] = (int8_t) _mm256_extract_epi8_rpl (vH, 31) - bias;
        }
#endif

        {
            __m256i vCompare = _mm256_cmpgt_epi8(vMaxH, vMaxHUnit);
            if (_mm256_movemask_epi8(vCompare)) {
                score = _mm256_hmax_epi8_rpl(vMaxH);
                /* if score has potential to overflow, abort early */
                if (score > maxp) {
                    result->flag |= PARASAIL_FLAG_SATURATED;
                    break;
                }
                vMaxHUnit = _mm256_set1_epi8(score);
                end_ref = j;
            }
        }

        /*if (score == stop) break;*/
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m256i vH = _mm256_load_si256(pvHStore+i);
        arr_store_col(result->rowcols->score_col, vH, i, segLen, bias);
    }
#endif

    if (score == INT8_MAX) {
        result->flag |= PARASAIL_FLAG_SATURATED;
    }

    if (parasail_result_is_saturated(result)) {
        score = INT8_MAX;
        end_query = 0;
        end_ref = 0;
    }
    else {
        if (end_ref == j-1) {
            /* end_ref was the last store column */
            __m256i *pv = pvHMax;
            pvHMax = pvHStore;
            pvHStore = pv;
        }
        else if (end_ref == j-2) {
            /* end_ref was the last load column */
            __m256i *pv = pvHMax;
            pvHMax = pvHLoad;
            pvHLoad = pv;
        }
        /* Trace the alignment ending position on read. */
        {
            int8_t *t = (int8_t*)pvHMax;
            int32_t column_len = segLen * segWidth;
            end_query = s1Len - 1;
            for (i = 0; i<column_len; ++i, ++t) {
                if (*t == score) {
                    int32_t temp = i / segWidth + i % segWidth * segLen;
                    if (temp < end_query) {
                        end_query = temp;
                    }
                }
            }
        }
    }

    result->score = score - bias;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(pvE);
    parasail_free(pvHMax);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

void parasail_profile_free(parasail_profile_t *profile)
{
    if (!profile) {
        fprintf(stderr, "%s: attempted free of NULL profile pointer\n", __func__);
        return;
    }

    if (NULL != profile->profile8.score) {
        profile->free(profile->profile8.score);
    }
    if (NULL != profile->profile8.matches) {
        profile->free(profile->profile8.matches);
    }
    if (NULL != profile->profile8.similar) {
        profile->free(profile->profile8.similar);
    }

    if (NULL != profile->profile16.score) {
        profile->free(profile->profile16.score);
    }
    if (NULL != profile->profile16.matches) {
        profile->free(profile->profile16.matches);
    }
    if (NULL != profile->profile16.similar) {
        profile->free(profile->profile16.similar);
    }

    if (NULL != profile->profile32.score) {
        profile->free(profile->profile32.score);
    }
    if (NULL != profile->profile32.matches) {
        profile->free(profile->profile32.matches);
    }
    if (NULL != profile->profile32.similar) {
        profile->free(profile->profile32.similar);
    }

    if (NULL != profile->profile64.score) {
        profile->free(profile->profile64.score);
    }
    if (NULL != profile->profile64.matches) {
        profile->free(profile->profile64.matches);
    }
    if (NULL != profile->profile64.similar) {
        profile->free(profile->profile64.similar);
    }

    free(profile);
}

void parasail_result_free(parasail_result_t *result)
{
    /* validate inputs */
    if (!result) {
        fprintf(stderr, "%s: attempted free of NULL result pointer\n", __func__);
        return;
    }
    
    if (result->flag & PARASAIL_FLAG_STATS) {
        if (result->flag & PARASAIL_FLAG_TABLE) {
            free(result->stats->tables->score_table);
            free(result->stats->tables->matches_table);
            free(result->stats->tables->similar_table);
            free(result->stats->tables->length_table);
            free(result->stats->tables);
        }
        if (result->flag & PARASAIL_FLAG_ROWCOL) {
            free(result->stats->rowcols->score_row);
            free(result->stats->rowcols->matches_row);
            free(result->stats->rowcols->similar_row);
            free(result->stats->rowcols->length_row);
            free(result->stats->rowcols->score_col);
            free(result->stats->rowcols->matches_col);
            free(result->stats->rowcols->similar_col);
            free(result->stats->rowcols->length_col);
            free(result->stats->rowcols);
        }
        free(result->stats);
    }
    else {
        if (result->flag & PARASAIL_FLAG_TABLE) {
            free(result->tables->score_table);
            free(result->tables);
        }
        if (result->flag & PARASAIL_FLAG_ROWCOL) {
            free(result->rowcols->score_row);
            free(result->rowcols->score_col);
            free(result->rowcols);
        }
    }
    if (result->flag & PARASAIL_FLAG_TRACE) {
        parasail_free(result->trace->trace_table);
        if (NULL != result->trace->trace_ins_table)
            parasail_free(result->trace->trace_ins_table);
        if (NULL != result->trace->trace_del_table)
            parasail_free(result->trace->trace_del_table);
        free(result->trace);
    }

    free(result);
}