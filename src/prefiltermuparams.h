#pragma once

// pattern=1101011

static const uint k = 5;
static const uint K = 7;
static const uint8_t s_Offsets[] = { 0, 1, 3, 5, 6 };
static const uint8_t *s_ptrOffsets = s_Offsets;
static const uint DICT_SIZE = 60466176;	// 36^5
static const uint RSB_SIZE = 500;
static const uint ALPHABET_SIZE = 36;

// N=60466176, Min=20, LoQ=43, Med=47, HiQ=51, Max=75, Avg=47.3611
static const int MIN_KMER_PAIR_SCORE = 25;

// Not used
static const int BIAS_WINDOW = 30; // 40
static const float TBIAS_SCALE = 0; // 0.15f;
static const float QBIAS_SCALE = 0; // 0.15f;
