#pragma once

static const char *prefiltermu_pattern = "1110011";
static const uint PREFILTER_KMER_NR_ONES = 5;
static const uint PREFILTER_KMER_WIDTH = 7;
static const uint8_t PREFILTER_KMER_ONES_OFFSETS[] = { 0, 1, 2, 5, 6 };
static const uint PREFILTER_KMER_DICT_SIZE = 60466176;	// 36^5
static const uint MAX_PREFILTER_KMER_HOOD_SIZE = 41293; // empirical cmd_kmrnbh()

// Can be overriden with -idxq or -idxt
static const uint MAX_QUERY_CHAINS_FOR_QUERY_NEIGHBORHOOD = 100;

extern bool g_QueryNeighborhood;