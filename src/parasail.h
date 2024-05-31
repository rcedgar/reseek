#pragma once

/***
This source code is partly derived from the parasail library
https://github.com/jeffdaily/parasail
Author jeffrey.daily@gmail.com
Copyright (c) 2015 Battelle Memorial Institute.
***/

#include <immintrin.h>
/*                                            3         2         1          */
/*                                           10987654321098765432109876543210*/
#define PARASAIL_FLAG_NW          (1 << 0) /*00000000000000000000000000000001*/
#define PARASAIL_FLAG_SG          (1 << 1) /*00000000000000000000000000000010*/
#define PARASAIL_FLAG_SW          (1 << 2) /*00000000000000000000000000000100*/
#define PARASAIL_FLAG_SG_S1_BEG   (1 << 3) /*00000000000000000000000000001000*/
#define PARASAIL_FLAG_SG_S1_END   (1 << 4) /*00000000000000000000000000010000*/
#define PARASAIL_FLAG_SATURATED   (1 << 6) /*00000000000000000000000001000000*/
#define PARASAIL_FLAG_BANDED      (1 << 7) /*00000000000000000000000010000000*/
#define PARASAIL_FLAG_NOVEC       (1 << 8) /*00000000000000000000000100000000*/
#define PARASAIL_FLAG_NOVEC_SCAN  (1 << 9) /*00000000000000000000001000000000*/
#define PARASAIL_FLAG_SCAN        (1 <<10) /*00000000000000000000010000000000*/
#define PARASAIL_FLAG_STRIPED     (1 <<11) /*00000000000000000000100000000000*/
#define PARASAIL_FLAG_DIAG        (1 <<12) /*00000000000000000001000000000000*/
#define PARASAIL_FLAG_BLOCKED     (1 <<13) /*00000000000000000010000000000000*/
#define PARASAIL_FLAG_SG_S2_BEG   (1 <<14) /*00000000000000000100000000000000*/
#define PARASAIL_FLAG_SG_S2_END   (1 <<15) /*00000000000000001000000000000000*/
#define PARASAIL_FLAG_STATS       (1 <<16) /*00000000000000010000000000000000*/
#define PARASAIL_FLAG_TABLE       (1 <<17) /*00000000000000100000000000000000*/
#define PARASAIL_FLAG_ROWCOL      (1 <<18) /*00000000000001000000000000000000*/
#define PARASAIL_FLAG_TRACE       (1 <<19) /*00000000000010000000000000000000*/
#define PARASAIL_FLAG_BITS_8      (1 <<20) /*00000000000100000000000000000000*/
#define PARASAIL_FLAG_BITS_16     (1 <<21) /*00000000001000000000000000000000*/
#define PARASAIL_FLAG_BITS_32     (1 <<22) /*00000000010000000000000000000000*/
#define PARASAIL_FLAG_BITS_64     (1 <<23) /*00000000100000000000000000000000*/
#define PARASAIL_FLAG_LANES_1     (1 <<24) /*00000001000000000000000000000000*/
#define PARASAIL_FLAG_LANES_2     (1 <<25) /*00000010000000000000000000000000*/
#define PARASAIL_FLAG_LANES_4     (1 <<26) /*00000100000000000000000000000000*/
#define PARASAIL_FLAG_LANES_8     (1 <<27) /*00001000000000000000000000000000*/
#define PARASAIL_FLAG_LANES_16    (1 <<28) /*00010000000000000000000000000000*/
#define PARASAIL_FLAG_LANES_32    (1 <<29) /*00100000000000000000000000000000*/
#define PARASAIL_FLAG_LANES_64    (1 <<30) /*01000000000000000000000000000000*/
#define PARASAIL_FLAG_INVALID  0x80000020  /*10000000000000000000000000100000*/

#define PARASAIL_MATRIX_TYPE_SQUARE 0
#define PARASAIL_MATRIX_TYPE_PSSM 1

#define restrict __restrict
#define PARASAIL_CHECK_NULL(x)	asserta((x) != 0)
#define PARASAIL_CHECK_GT0(x)	asserta((x) > 0)
#define PARASAIL_CHECK_GE0(x)	asserta((x) >= 0)

#define SWAP(A,B) { __m256i* tmp = A; A = B; B = tmp; }
#define SWAP3(A,B,C) { __m256i* tmp = A; A = B; B = C; C = tmp; }
#define NEG_INF INT8_MIN

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)
#define _mm256_insert_epi8_rpl _mm256_insert_epi8

/* for traceback */
#define PARASAIL_ZERO_MASK 120 /* all bits set except the first three */
#define PARASAIL_E_MASK 103 /* all bits set except the E bits */
#define PARASAIL_F_MASK 31 /* all bits set except the F bits */
#define PARASAIL_ZERO   0
#define PARASAIL_INS    1
#define PARASAIL_DEL    2
#define PARASAIL_DIAG   4
#define PARASAIL_DIAG_E 8
#define PARASAIL_INS_E  16
#define PARASAIL_DIAG_F 32
#define PARASAIL_DEL_F  64

typedef struct parasail_matrix {
    const char * name;
    const int * matrix;
    const int * mapper;
    int size;
    int max;
    int min;
    int * user_matrix;
    int type;
    int length;
    const char * alphabet;
    const char * query;
} parasail_matrix_t;

typedef struct parasail_profile_data {
    void * score;
    void * matches;
    void * similar;
} parasail_profile_data_t;

typedef struct parasail_profile {
    const char *s1;
    int s1Len;
    const parasail_matrix_t *matrix;
    struct parasail_profile_data profile8;
    struct parasail_profile_data profile16;
    struct parasail_profile_data profile32;
    struct parasail_profile_data profile64;
    void (*free)(void * profile);
    int stop;
} parasail_profile_t;

typedef union __m256i_8 {
    __m256i m;
    int8_t v[32];
} __m256i_8_t;


typedef struct parasail_result_extra_stats_tables {
    int * restrict score_table;     /* DP table of scores */
    int * restrict matches_table;   /* DP table of exact match counts */
    int * restrict similar_table;   /* DP table of similar substitution counts */
    int * restrict length_table;    /* DP table of lengths */
} parasail_result_extra_stats_tables_t;

typedef struct parasail_result_extra_stats_rowcols {
    int * restrict score_row;       /* last row of DP table of scores */
    int * restrict matches_row;     /* last row of DP table of exact match counts */
    int * restrict similar_row;     /* last row of DP table of similar substitution counts */
    int * restrict length_row;      /* last row of DP table of lengths */
    int * restrict score_col;       /* last col of DP table of scores */
    int * restrict matches_col;     /* last col of DP table of exact match counts */
    int * restrict similar_col;     /* last col of DP table of similar substitution counts */
    int * restrict length_col;      /* last col of DP table of lengths */
} parasail_result_extra_stats_rowcols_t;

typedef struct parasail_result_extra_stats {
    int matches;    /* number of exactly matching characters in the alignment */
    int similar;    /* number of similar characters (positive substitutions) in the alignment */
    int length;     /* length of the alignment */
    union {
        void *extra;
        parasail_result_extra_stats_tables_t *tables;
        parasail_result_extra_stats_rowcols_t *rowcols;
    };
} parasail_result_extra_stats_t;

typedef struct parasail_result_extra_tables {
    int * restrict score_table;     /* DP table of scores */
} parasail_result_extra_tables_t;

typedef struct parasail_result_extra_rowcols {
    int * restrict score_row;       /* last row of DP table of scores */
    int * restrict score_col;       /* last col of DP table of scores */
} parasail_result_extra_rowcols_t;

typedef struct parasail_result_extra_trace {
    void * restrict trace_table;    /* DP table of traceback */
    void * restrict trace_ins_table;/* DP table of insertions traceback */
    void * restrict trace_del_table;/* DP table of deletions traceback */
} parasail_result_extra_trace_t;

typedef struct parasail_result {
    int score;      /* alignment score */
    int end_query;  /* end position of query sequence */
    int end_ref;    /* end position of reference sequence */
    int flag;       /* bit field for various flags */
    /* union of pointers to extra result data based on the flag */
    union {
        void *extra;
        parasail_result_extra_stats_t *stats;
        parasail_result_extra_tables_t *tables;
        parasail_result_extra_rowcols_t *rowcols;
        parasail_result_extra_trace_t *trace;
    };
} parasail_result_t;

typedef struct parasail_cigar_ {
    uint32_t *seq;
    int len;
    int64_t beg_query;
    int64_t beg_ref;
} parasail_cigar_t;

static inline __m256i arr_load(
        __m256i *array,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    return _mm256_load_si256(array + (1LL*d*seglen+t));
}

static inline __m256i arr_load(
        void *array,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    return _mm256_load_si256((__m256i *) array + (1LL*d*seglen+t));
}

static inline void arr_store(
        __m256i *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    _mm256_store_si256(array + (1LL*d*seglen+t), vH);
}

static inline void arr_store(
        void *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    _mm256_store_si256((__m256i *) array + (1LL*d*seglen+t), vH);
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

void* parasail_memalign(size_t alignment, size_t size);


static inline int8_t _mm256_hmax_epi8_rpl(__m256i a) {
    a = _mm256_max_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 8));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 4));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 2));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 1));
    return _mm256_extract_epi8(a, 31);
}

static inline void parasail_memset___m256i(__m256i *b, __m256i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm256_store_si256(&b[i], c);
    }
}

void parasail_free(void *ptr);

static inline uint32_t* parasail_reverse_uint32_t(const uint32_t *s, size_t length)
{
    uint32_t *r = NULL;
    size_t i = 0;
    size_t j = 0;

    PARASAIL_CALLOC(r, uint32_t, length);
    for (i=0,j=length-1; i<length; ++i,--j) {
        r[i] = s[j];
    }

    return r;
}

static inline __m256i * parasail_memalign___m256i(size_t alignment, size_t size)
{
    return (__m256i *) parasail_memalign(alignment, size*sizeof(__m256i));
}

static inline void parasail_free___m256i(void *ptr)
{
    parasail_free((__m256i*)ptr);
}

static inline int parasail_result_is_trace(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_TRACE;
}

char* parasail_cigar_decode(parasail_cigar_t *cigar);
parasail_cigar_t* parasail_result_get_cigar_extra(
        parasail_result_t *result,
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const parasail_matrix_t *matrix,
        int case_sensitive,
        const char *alphabet_aliases);

/////////////////////////////////////////////////////////////////////////////
parasail_result_t* parasail_result_new();

parasail_result_t* parasail_result_new_trace(const int a, const int b,
  const size_t alignment, const size_t size);

parasail_profile_t * parasail_profile_create_avx_256_8(
        const char * const restrict s1, const int _s1Len,
        const parasail_matrix_t *matrix);

void parasail_profile_free(parasail_profile_t *profile);

parasail_result_t* parasail_sw_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

parasail_result_t* parasail_sw_trace_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);
