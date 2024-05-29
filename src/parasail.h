#pragma once

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
