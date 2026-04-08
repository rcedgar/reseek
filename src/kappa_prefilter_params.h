#pragma once

static const char *prefilter_kappa_pattern = "1110011";
static const uint KAPPA_PREFILTER_KMER_NR_ONES = 5;
static const uint KAPPA_PREFILTER_KMER_WIDTH = 7;
static const uint8_t KAPPA_PREFILTER_KMER_ONES_OFFSETS[] = { 0, 1, 2, 5, 6 };
//static const uint KAPPA_PREFILTER_KMER_DICT_SIZE = 60466176;	// 36^5 (for Mu)
static const uint KAPPA_PREFILTER_KMER_DICT_SIZE = 33554432;	// 32^5
static const uint MAX_KAPPA_PREFILTER_KMER_HOOD_SIZE = 83940; // empirical cmd_kappa_kmrnbh()
static const uint KAPPA_AS = 32;