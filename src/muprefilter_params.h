#pragma once

// pattern=1110011
static const char *muprefilter_pattern = "1110011";
static const uint k = 5;
static const uint K = 7;
static const uint8_t s_Offsets[] = { 0, 1, 2, 5, 6 };

static const uint8_t *s_ptrOffsets = s_Offsets;
static const uint DICT_SIZE = 60466176;	// 36^5
static const uint RSB_SIZE = 1500;
static const uint ALPHABET_SIZE = 36;

// N=1679616, Min=16, LoQ=34, Med=38, HiQ=41, Max=60, Avg=37.8889
//static const int MIN_KMER_PAIR_SCORE = 36;
static const int MIN_KMER_PAIR_SCORE = 29;
