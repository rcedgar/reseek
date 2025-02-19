#pragma once

#define	USE_BIAS		0
#define	KMER_SORT		0
#define KMER_COMPLEXITY	0

// pattern=1110011
static const uint k = 5;
static const uint K = 7;
static const uint8_t s_Offsets[] = { 0, 1, 2, 5, 6 };

static const uint8_t *s_ptrOffsets = s_Offsets;
static const uint DICT_SIZE = 60466176;	// 36^5
static const uint RSB_SIZE = 1500;
static const uint ALPHABET_SIZE = 36;

// N=1679616, Min=16, LoQ=34, Med=38, HiQ=41, Max=60, Avg=37.8889
#if KMER_SORT
static const int MIN_KMER_PAIR_SCORE = 0;
#define QUERY_COUNT_MULTIPLIER	5000	// something wrong here, too slow
#else
static const int MIN_KMER_PAIR_SCORE = 36;
#endif

#if USE_BIAS
static const int BIAS_WINDOW = 30; // 40
static const float TBIAS_SCALE = 0; // 0.15f;
static const float QBIAS_SCALE = 0; // 0.15f;
#endif
