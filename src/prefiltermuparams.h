#pragma once

// pattern=1110011
static const uint k = 5;
static const uint K = 7;
static const uint8_t s_Offsets[] = { 0, 1, 2, 5, 6 };

static const uint8_t *s_ptrOffsets = s_Offsets;
static const uint DICT_SIZE = 60466176;	// 36^5
static const uint RSB_SIZE = 1500;
static const uint ALPHABET_SIZE = 36;

static const int MIN_KMER_PAIR_SCORE = 36;

// Not used
static const int BIAS_WINDOW = 30; // 40
static const float TBIAS_SCALE = 0; // 0.15f;
static const float QBIAS_SCALE = 0; // 0.15f;
