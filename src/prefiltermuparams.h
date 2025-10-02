#pragma once

#define	KMER_SORT		0
#define KMER_COMPLEXITY	0

// pattern=1110011
static const char *prefiltermu_pattern = "1110011";
static const uint k = 5;
static const uint K = 7;
static const uint8_t s_Offsets[] = { 0, 1, 2, 5, 6 };

static const uint8_t *s_ptrOffsets = s_Offsets;
static const uint DICT_SIZE = 60466176;	// 36^5
static const uint MAX_HOOD_SIZE = 41293; // empirical cmd_kmrnbh()
static const uint RSB_SIZE = 1500;
static const uint ALPHABET_SIZE = 36;
static const uint MAX_QUERY_CHAINS_FOR_QUERY_NEIGHBORHOOD = 100;

static const int MIN_KMER_PAIR_SCORE = 36;

extern bool g_QueryNeighborhood;