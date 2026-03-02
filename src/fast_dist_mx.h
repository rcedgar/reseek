#pragma once

#include <immintrin.h>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>

// max |i-j|
static const uint32_t dist_mx_band_width = 100;

// band_w(M) = 99
static const uint32_t band_w_M = dist_mx_band_width-1;

 // band_c(band_w(M)) = 4950
static const uint32_t band_c_w_M = (band_w_M * (band_w_M + 1)) / 2;

// base(i) = sum_{r=0}^{i-1} min(r,w)
static inline uint32_t band_base(uint32_t i) {
	if (i <= dist_mx_band_width) return (i * (i - 1)) / 2;
	return band_w_M * i - band_c_w_M;
}

static inline uint32_t band_K(uint32_t L) {
	return band_base(L);
}

static inline uint32_t band_ij_to_k(uint32_t i, uint32_t j, uint32_t L) {
	(void)L;
	uint32_t w = band_w_M;
	assert(j < i);
	assert(i - j <= w);

	uint32_t len   = std::min(i, w);
	uint32_t start = i - len;
	return band_base(i) + (j - start);
}

static inline void banded_distances_avx2_u16_xyz(
	const uint16_t* __restrict xyz,   // length 3*L: [x...][y...][z...], each uint16 IC
	uint32_t L,
	uint16_t* __restrict outIC        // length band_K(L, bi)
){
	asserta(xyz && outIC);
	asserta(L >= 8);

	const uint16_t* __restrict x = xyz;
	const uint16_t* __restrict y = xyz + L;
	const uint16_t* __restrict z = xyz + 2u*L;

	const __m256i i10000 = _mm256_set1_epi32(10000);
	const __m256i izero  = _mm256_setzero_si256();
	const __m256i i65535 = _mm256_set1_epi32(65535);

	for (uint32_t i = 1; i < L; ++i) {
		uint32_t len = std::min(i, band_w_M);
		if (len == 0) continue;

		uint32_t start = i - len;		// max(0, i-w)
		uint32_t base  = band_base(i);	// row base in flattened output
		uint32_t t = 0;

		const __m256i xi = _mm256_set1_epi32((int32_t)x[i]);
		const __m256i yi = _mm256_set1_epi32((int32_t)y[i]);
		const __m256i zi = _mm256_set1_epi32((int32_t)z[i]);

		for (; t + 7 < len; t += 8) {
			uint32_t j = start + t;

			__m128i x8u16 = _mm_loadu_si128((const __m128i*)(x + j));
			__m128i y8u16 = _mm_loadu_si128((const __m128i*)(y + j));
			__m128i z8u16 = _mm_loadu_si128((const __m128i*)(z + j));

			__m256i xj = _mm256_cvtepu16_epi32(x8u16);
			__m256i yj = _mm256_cvtepu16_epi32(y8u16);
			__m256i zj = _mm256_cvtepu16_epi32(z8u16);

			__m256i dx = _mm256_sub_epi32(xj, xi);
			__m256i dy = _mm256_sub_epi32(yj, yi);
			__m256i dz = _mm256_sub_epi32(zj, zi);

			__m256i d2 = _mm256_mullo_epi32(dx, dx);
			d2 = _mm256_add_epi32(d2, _mm256_mullo_epi32(dy, dy));
			d2 = _mm256_add_epi32(d2, _mm256_mullo_epi32(dz, dz));

			__m256 fs = _mm256_sqrt_ps(_mm256_cvtepi32_ps(d2));

			// deterministic round-to-nearest: round then truncate
			fs = _mm256_round_ps(fs, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
			__m256i s_int = _mm256_cvttps_epi32(fs);

			__m256i dic = _mm256_add_epi32(s_int, i10000);
			dic = _mm256_max_epi32(dic, izero);
			dic = _mm256_min_epi32(dic, i65535);

			// pack 8x i32 -> 8x u16 contiguously (correct)
			__m128i dic_lo = _mm256_castsi256_si128(dic);
			__m128i dic_hi = _mm256_extracti128_si256(dic, 1);
			__m128i out8u16 = _mm_packus_epi32(dic_lo, dic_hi);

			_mm_storeu_si128((__m128i*)(outIC + base + t), out8u16);
		}

		for (; t < len; ++t) {
			uint32_t j = start + t;

			int32_t dx = (int32_t)x[j] - (int32_t)x[i];
			int32_t dy = (int32_t)y[j] - (int32_t)y[i];
			int32_t dz = (int32_t)z[j] - (int32_t)z[i];
			int32_t d2 = dx*dx + dy*dy + dz*dz;

			float s = std::sqrt((float)d2);                 // deci-units
			uint32_t si = (uint32_t)std::floor(s + 0.5f);   // round-to-nearest (s>=0)
			uint32_t dic = si + 10000u;
			if (dic > 65535u) dic = 65535u;
			outIC[base + t] = (uint16_t)dic;
		}
	}
}
