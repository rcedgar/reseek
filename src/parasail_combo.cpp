#include "myutils.h"
#include "parasail.h"
#include "dssaligner.h"
#include "timing.h"

#define restrict __restrict
#define PARASAIL_CHECK_NULL(x)	asserta((x) != 0)
#define PARASAIL_CHECK_GT0(x)	asserta((x) > 0)
#define PARASAIL_CHECK_GE0(x)	asserta((x) >= 0)

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)

//static const int parasail_combo_map[256] = {
//0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
//15,16,17,18,19,20,21,22,0,24,25,0,0,0,0,0,0,26,27,28,29,30,31,32,33,34,35,0,0,0,0,0,0,0,0,0,0,0,0,0,
//0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
//0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};

static const int parasail_combo_map[256] = {
0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};

static const int parasail_combo_[36*36] = {
 3,-1, 1, 1,-3,-1, 2,-2, 0, 1,-3,-1,-1,-5,-3, 0,-4,-2, 0,-4,-2,-2,-6,-4,-1,-5,-3,-1,-5,-3,-3,-7,-5,-2,-6,-4,  // 0
-1, 4, 2,-3, 1,-1,-2, 2, 0,-3, 2, 0,-5, 0,-2,-4, 0,-2,-4, 1,-1,-6,-1,-3,-5,-1,-3,-5, 0,-2,-7,-2,-4,-6,-2,-4,  // 1
 1, 2, 3,-1,-1, 0, 0, 0, 1,-1, 0, 1,-3,-2,-1,-2,-2,-1,-2,-1, 0,-4,-3,-2,-3,-3,-2,-3,-2,-1,-5,-4,-3,-4,-4,-3,  // 2
 1,-3,-1, 3, 0, 1, 2,-2, 0,-1,-5,-3, 2,-2,-1, 0,-4,-2,-2,-6,-4, 0,-3,-2,-1,-5,-3,-3,-7,-5,-1,-4,-3,-2,-6,-4,  // 3
-3, 1,-1, 0, 4, 2,-2, 3, 1,-5, 0,-2,-2, 2, 0,-4, 1,-1,-6,-1,-3,-3, 1,-1,-5, 0,-2,-7,-2,-4,-4, 0,-2,-6,-1,-3,  // 4
-1,-1, 0, 1, 2, 3, 0, 1, 2,-3,-2,-1,-1, 0, 1,-2,-1, 0,-4,-3,-2,-2,-1, 0,-3,-2,-1,-5,-4,-3,-3,-2,-1,-4,-3,-2,  // 5
 2,-2, 0, 2,-2, 0, 3,-1, 1, 0,-4,-2, 0,-4,-2, 1,-3,-1,-1,-5,-3,-1,-5,-3, 0,-4,-2,-2,-6,-4,-2,-6,-4,-1,-5,-3,  // 6
-2, 2, 0,-2, 3, 1,-1, 3, 1,-4, 0,-2,-4, 1,-1,-3, 1,-1,-5,-1,-3,-5, 0,-2,-4, 0,-2,-6,-2,-4,-6,-1,-3,-5,-1,-3,  // 7
 0, 0, 1, 0, 1, 2, 1, 1, 2,-2,-2,-1,-2,-1, 0,-1,-1, 0,-3,-3,-2,-3,-2,-1,-2,-2,-1,-4,-4,-3,-4,-3,-2,-3,-3,-2,  // 8
 1,-3,-1,-1,-5,-3, 0,-4,-2, 2,-2, 0, 0,-4,-2, 1,-3,-1, 1,-3,-1,-1,-5,-3, 0,-4,-2, 1,-3,-1,-2,-5,-4,-1,-5,-3,  // 9
-3, 2, 0,-5, 0,-2,-4, 0,-2,-2, 3, 1,-4, 0,-2,-3, 1,-1,-3, 2, 0,-5, 0,-2,-4, 0,-2,-3, 1,-1,-5,-1,-3,-5, 0,-2,  // 10
-1, 0, 1,-3,-2,-1,-2,-2,-1, 0, 1, 2,-2,-2,-1,-1,-1, 0,-1, 0, 1,-3,-2,-1,-2,-2,-1,-1,-1, 0,-4,-3,-2,-3,-2,-1,  // 11
-1,-5,-3, 2,-2,-1, 0,-4,-2, 0,-4,-2, 2,-1, 0, 1,-3,-1,-1,-5,-3, 1,-2,-1, 0,-4,-2,-2,-5,-4, 1,-3,-1, 0,-4,-2,  // 12
-5, 0,-2,-2, 2, 0,-4, 1,-1,-4, 0,-2,-1, 3, 1,-3, 2, 0,-5, 0,-2,-2, 2, 0,-4, 1,-1,-5,-1,-3,-3, 1,-1,-4, 0,-2,  // 13
-3,-2,-1,-1, 0, 1,-2,-1, 0,-2,-2,-1, 0, 1, 2,-1, 0, 1,-3,-2,-1,-1, 0, 1,-2,-1, 0,-4,-3,-2,-1,-1, 0,-2,-2,-1,  // 14
 0,-4,-2, 0,-4,-2, 1,-3,-1, 1,-3,-1, 1,-3,-1, 2,-2, 0, 0,-4,-2, 0,-4,-2, 1,-3,-1,-1,-5,-3, 0,-4,-2, 0,-3,-2,  // 15
-4, 0,-2,-4, 1,-1,-3, 1,-1,-3, 1,-1,-3, 2, 0,-2, 2, 0,-4, 0,-2,-4, 1,-1,-3, 1,-1,-5, 0,-2,-4, 0,-2,-3, 1,-1,  // 16
-2,-2,-1,-2,-1, 0,-1,-1, 0,-1,-1, 0,-1, 0, 1, 0, 0, 1,-2,-2,-1,-2,-1, 0,-1,-1, 0,-3,-2,-1,-2,-2,-1,-2,-1, 0,  // 17
 0,-4,-2,-2,-6,-4,-1,-5,-3, 1,-3,-1,-1,-5,-3, 0,-4,-2, 3,-1, 1, 0,-3,-2, 1,-2,-1, 1,-2,-1,-1,-5,-3, 0,-4,-2,  // 18
-4, 1,-1,-6,-1,-3,-5,-1,-3,-3, 2, 0,-5, 0,-2,-4, 0,-2,-1, 3, 1,-3, 1,-1,-2, 2, 0,-2, 2, 0,-5, 0,-2,-4, 0,-2,  // 19
-2,-1, 0,-4,-3,-2,-3,-3,-2,-1, 0, 1,-3,-2,-1,-2,-2,-1, 1, 1, 2,-2,-1, 0,-1, 0, 1,-1, 0, 1,-3,-2,-1,-2,-2,-1,  // 20
-2,-6,-4, 0,-3,-2,-1,-5,-3,-1,-5,-3, 1,-2,-1, 0,-4,-2, 0,-3,-2, 3,-1, 1, 2,-2, 0,-1,-5,-3, 2,-2, 0, 0,-4,-2,  // 21
-6,-1,-3,-3, 1,-1,-5, 0,-2,-5, 0,-2,-2, 2, 0,-4, 1,-1,-3, 1,-1,-1, 4, 2,-2, 2, 0,-5, 0,-2,-2, 2, 0,-4, 1,-1,  // 22
-4,-3,-2,-2,-1, 0,-3,-2,-1,-3,-2,-1,-1, 0, 1,-2,-1, 0,-2,-1, 0, 1, 2, 3, 0, 0, 1,-3,-2,-1, 0, 0, 1,-2,-1, 0,  // 23
-1,-5,-3,-1,-5,-3, 0,-4,-2, 0,-4,-2, 0,-4,-2, 1,-3,-1, 1,-2,-1, 2,-2, 0, 2,-1, 0, 0,-4,-2, 0,-4,-2, 1,-3,-1,  // 24
-5,-1,-3,-5, 0,-2,-4, 0,-2,-4, 0,-2,-4, 1,-1,-3, 1,-1,-2, 2, 0,-2, 2, 0,-1, 3, 1,-4, 0,-2,-4, 1,-1,-3, 2, 0,  // 25
-3,-3,-2,-3,-2,-1,-2,-2,-1,-2,-2,-1,-2,-1, 0,-1,-1, 0,-1, 0, 1, 0, 0, 1, 0, 1, 2,-2,-2,-1,-2,-1, 0,-1, 0, 1,  // 26
-1,-5,-3,-3,-7,-5,-2,-6,-4, 1,-3,-1,-2,-5,-4,-1,-5,-3, 1,-2,-1,-1,-5,-3, 0,-4,-2, 3,-1, 1, 0,-3,-2, 1,-3,-1,  // 27
-5, 0,-2,-7,-2,-4,-6,-2,-4,-3, 1,-1,-5,-1,-3,-5, 0,-2,-2, 2, 0,-5, 0,-2,-4, 0,-2,-1, 3, 1,-3, 1,-1,-3, 2, 0,  // 28
-3,-2,-1,-5,-4,-3,-4,-4,-3,-1,-1, 0,-4,-3,-2,-3,-2,-1,-1, 0, 1,-3,-2,-1,-2,-2,-1, 1, 1, 2,-2,-1, 0,-1, 0, 1,  // 29
-3,-7,-5,-1,-4,-3,-2,-6,-4,-2,-5,-4, 1,-3,-1, 0,-4,-2,-1,-5,-3, 2,-2, 0, 0,-4,-2, 0,-3,-2, 3,-1, 1, 2,-2,-1,  // 30
-7,-2,-4,-4, 0,-2,-6,-1,-3,-5,-1,-3,-3, 1,-1,-4, 0,-2,-5, 0,-2,-2, 2, 0,-4, 1,-1,-3, 1,-1,-1, 3, 1,-2, 2, 0,  // 31
-5,-4,-3,-3,-2,-1,-4,-3,-2,-4,-3,-2,-1,-1, 0,-2,-2,-1,-3,-2,-1, 0, 0, 1,-2,-1, 0,-2,-1, 0, 1, 1, 2,-1, 0, 1,  // 32
-2,-6,-4,-2,-6,-4,-1,-5,-3,-1,-5,-3, 0,-4,-2, 0,-3,-2, 0,-4,-2, 0,-4,-2, 1,-3,-1, 1,-3,-1, 2,-2,-1, 2,-2, 0,  // 33
-6,-2,-4,-6,-1,-3,-5,-1,-3,-5, 0,-2,-4, 0,-2,-3, 1,-1,-4, 0,-2,-4, 1,-1,-3, 2, 0,-3, 2, 0,-2, 2, 0,-2, 3, 1,  // 34
-4,-4,-3,-4,-3,-2,-3,-3,-2,-3,-2,-1,-2,-2,-1,-2,-1, 0,-2,-2,-1,-2,-1, 0,-1, 0, 1,-1, 0, 1,-1, 0, 1, 0, 1, 2,  // 35
};

parasail_matrix_t parasail_combo_matrix = {
	"combo",  // name
	parasail_combo_,  // matrix
	parasail_combo_map,  // mapper
	36,  // size
	4, // max
	-7,  // min
	NULL, // user_matrix
	PARASAIL_MATRIX_TYPE_SQUARE,  // type
	36,  // length
	"ABCDEFGHIJKLMNOPQRSTUVWZYZabcdefghij", // alphabet
	NULL // query
};

static inline int8_t _mm256_hmax_epi8_rpl(__m256i a) {
    a = _mm256_max_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 8));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 4));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 2));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 1));
    return _mm256_extract_epi8(a, 31);
}

static void parasail_memset___m256i(__m256i *b, __m256i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm256_store_si256(&b[i], c);
    }
}

#define PARASAIL_CALLOC(var,type,count) do {                                    \
    size_t _size = sizeof(type)*(count);                                        \
    var = (type*)malloc(_size);                                                 \
    if (!var) {                                                                 \
        fprintf(stderr, "%s: failed to malloc %zu bytes\n", __func__, (_size)); \
        return NULL;                                                            \
    }                                                                           \
} while(0)

#define PARASAIL_NEW(var,type) PARASAIL_CALLOC(var,type,1)

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
    ptr = _aligned_malloc(size, alignment);
	asserta(ptr != 0);
    return ptr;
}

static __m256i * parasail_memalign___m256i(size_t alignment, size_t size)
{
    return (__m256i *) parasail_memalign(alignment, size*sizeof(__m256i));
}

void parasail_free(void *ptr)
{
     _aligned_free(ptr);
}

static void parasail_free___m256i(void *ptr)
{
    parasail_free((__m256i*)ptr);
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

typedef struct parasail_result {
    int score;      /* alignment score */
    int end_query;  /* end position of query sequence */
    int end_ref;    /* end position of reference sequence */
    int flag;       /* bit field for various flags */
} parasail_result_t;

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

static int parasail_result_is_saturated(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_SATURATED;
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

float DSSAligner::AlignCombo_Prof_Para(const vector<byte> &LettersA,
  const vector<byte> &LettersB)
	{
	uint LA = SIZE(LettersA);
	uint LB = SIZE(LettersB);
	m_ComboLettersA = &LettersA;
	m_ComboLettersB = &LettersB;
	//SetProf_Combo_Para();
	StartTimer(SWPara);

	const parasail_profile_t * const restrict profile =
	  (const parasail_profile_t * const restrict) m_ProfPara;
	parasail_result_t* result =
	  parasail_sw_striped_profile_avx2_256_8(profile, (const char *) LettersB.data(), LB, 99, 99);

	const parasail_profile_t * const restrict profile_rev =
	  (const parasail_profile_t * const restrict) m_ProfParaRev;
	parasail_result_t* result_rev =
	  parasail_sw_striped_profile_avx2_256_8(profile_rev, (const char *) LettersB.data(), LB, 99, 99);

	EndTimer(SWPara);
	if (result->flag & PARASAIL_FLAG_SATURATED)
		result->score = 777;
	if (result_rev->flag & PARASAIL_FLAG_SATURATED)
		result_rev->score = 777;
	float Score = (float) result->score - (float) result_rev->score;
	free(result);
	free(result_rev);
	return Score;
	}

void DSSAligner::SetProf_Combo_Para()
	{
	StartTimer(SetProf_Combo_Para);
	if (m_ProfPara != 0)
		free(m_ProfPara);
	if (m_ProfParaRev != 0)
		free(m_ProfParaRev);
	const char * A = (const char *) m_ComboLettersA->data();
	uint LA = SIZE(*m_ComboLettersA);
	m_ProfPara = parasail_profile_create_avx_256_8(A, LA, &parasail_combo_matrix);

	vector<byte> ARev = *m_ComboLettersA;
	reverse(ARev.begin(), ARev.end());
	m_ProfParaRev = parasail_profile_create_avx_256_8((const char *) ARev.data(), LA, &parasail_combo_matrix);
	EndTimer(SetProf_Combo_Para);
	}
