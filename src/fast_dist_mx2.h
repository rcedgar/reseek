#pragma once

#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>

struct BandIndexLite {
    uint32_t M = 0;      // abs(i-j) < M
    uint32_t w = 0;      // w = M-1, max stored offset in indices
    uint32_t c = 0;      // c = w*(w+1)/2

    explicit BandIndexLite(uint32_t M_) : M(M_), w(M_ > 0 ? M_ - 1 : 0), c(w * (w + 1) / 2) {}
};

///////////////////////////////////////////////////////////////
// START alternative to BandIndexLite
static inline uint32_t band_w(uint32_t M) { return (M > 0 ? M - 1 : 0); }
static inline uint32_t band_c(uint32_t w) { return (w * (w + 1)) / 2; }

// base(i) = sum_{r=0}^{i-1} min(r,w)
static inline uint32_t band_base(uint32_t i, uint32_t M) {
    uint32_t w = band_w(M);
    uint32_t c = band_c(w);
    if (i <= w + 1) return (i * (i - 1)) / 2;
    return w * i - c;
}

static inline uint32_t band_K(uint32_t N, uint32_t M) {
    return band_base(N, M);
}

static inline uint32_t band_ij_to_k(uint32_t i, uint32_t j, uint32_t N, uint32_t M) {
    (void)N;
    uint32_t w = band_w(M);
    assert(j < i);
    assert(i - j <= w);

    uint32_t len   = std::min(i, w);
    uint32_t start = i - len;
    return band_base(i, M) + (j - start);
}
// END alternative to BandIndexLite
///////////////////////////////////////////////////////////////

// base(i) = number of stored entries in rows 0..i-1
static inline uint32_t band_base(uint32_t i, const BandIndexLite& bi) {
    // sum_{r=0}^{i-1} min(r,w)
    // if i <= w+1: 0+1+...+(i-1) = i(i-1)/2
    // else:        c + (i-1-w)*w = w*i - c
    if (i <= bi.w + 1) return (i * (i - 1)) / 2;
    return bi.w * i - bi.c;
}

static inline uint32_t band_K(uint32_t N, const BandIndexLite& bi) {
    return band_base(N, bi);
}

static inline uint32_t band_ij_to_k(uint32_t i, uint32_t j, uint32_t N, const BandIndexLite& bi) {
    (void)N;
    assert(j < i);
    assert(i - j <= bi.w);

    uint32_t len   = std::min(i, bi.w);
    uint32_t start = i - len;              // max(0, i-w)
    return band_base(i, bi) + (j - start);
}

static inline void band_k_to_ij(uint32_t k, uint32_t N, const BandIndexLite& bi, uint32_t& i, uint32_t& j) {
    const uint32_t K = band_K(N, bi);
    assert(k < K);

    uint32_t ii;
    if (bi.w == 0) {
        // M=1 => no entries at all; caller shouldn't call
        ii = 0;
    } else if (k < bi.c) {
        // Solve ii(ii-1)/2 <= k < (ii+1)ii/2
        // ii ≈ floor((1 + sqrt(1+8k))/2)
        double r = std::sqrt(1.0 + 8.0 * (double)k);
        ii = (uint32_t)((1.0 + r) * 0.5);

        // fixup
        while (band_base(ii, bi) > k) --ii;
        while (ii + 1 <= N && band_base(ii + 1, bi) <= k) ++ii;
    } else {
        // base(ii) = w*ii - c  => ii ≈ floor((k + c)/w)
        ii = (k + bi.c) / bi.w;

        // fixup (rare; protects boundaries)
        while (band_base(ii, bi) > k) --ii;
        while (ii + 1 <= N && band_base(ii + 1, bi) <= k) ++ii;
    }

    i = ii;

    uint32_t len   = std::min(i, bi.w);
    uint32_t start = i - len;
    uint32_t t = k - band_base(i, bi);
    j = start + t;

    assert(j < i);
    assert(i - j <= bi.w);
}

static inline void banded_distances_avx2_u16_xyz(
	const uint16_t* __restrict xyz,   // length 3*N: [x...][y...][z...], each uint16 IC
	uint32_t N,
	const BandIndexLite& bi,          // constructed once for your fixed M (e.g. M=100)
	uint16_t* __restrict outIC        // length band_K(N, bi)
){
	assert(xyz && outIC);
	if (N < 2 || bi.w == 0) return;

	const uint16_t* __restrict x = xyz;
	const uint16_t* __restrict y = xyz + N;
	const uint16_t* __restrict z = xyz + 2u*N;

	const __m256i i10000 = _mm256_set1_epi32(10000);
	const __m256i izero  = _mm256_setzero_si256();
	const __m256i i65535 = _mm256_set1_epi32(65535);

	for (uint32_t i = 1; i < N; ++i) {
		uint32_t len = std::min(i, bi.w);
		if (len == 0) continue;

		uint32_t start = i - len;                 // max(0, i-w)
		uint32_t base  = band_base(i, bi);        // row base in flattened output
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
