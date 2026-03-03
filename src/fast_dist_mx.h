#pragma once

#include <immintrin.h>
#include "dss.h"

/***
Fixed-stride symmetric band for repeated random access.
Store a full band around the diagonal (both sides),
with a fixed stride per row. Then (i,j) -> k is just
a couple ops.
***/

static uint16_t CoordToIC(float X) { return uint16_t((X + 1000)*10 + 0.5); }
static float ICToCoord(uint16_t IC) { return float(IC/10.0f) - 1000; }

static uint16_t coord2ic(float X) { return uint16_t((X + 1000)*10 + 0.5); }
static float ic2coord(uint16_t IC) { return float(IC/10.0f) - 1000; }

static const uint32_t dist_mx_band_width = 100;
static const uint32_t nn_min_offset = 12;

// Let w = M-1 and stride = 2*w+1 (includes diagonal).

static const auto M = dist_mx_band_width;
static const auto w = M - 1;
static const uint32_t stride = 2*w + 1;

static inline uint32_t band_ij_to_k(uint32_t i, uint32_t j) {
	const uint32_t w = (M > 0 ? M - 1 : 0);
	const int32_t  d = (int32_t)j - (int32_t)i;
	// abs(i-j) < M  <=>  -w <= d <= w
	assert(d >= -(int32_t)w && d <= (int32_t)w);
	const uint32_t stride = 2*w + 1;
	return i*stride + (uint32_t)(d + (int32_t)w);
}

static inline void bandk_ij(uint32_t k, uint32_t& i, uint32_t& j) {
	i = k / stride;
	uint32_t off = k - i*stride;
	int32_t d = (int32_t)off - (int32_t)w;
	j = (uint32_t)((int32_t)i + d);
}

// dist mx flattened length K = L*(2*(M-1)+1)
static inline uint32_t band_K(uint32_t L) {
	return L*(2*(M-1)+1);
}

static inline uint16_t band2_get_checked(const uint16_t* outIC, uint32_t i, uint32_t j) {
    int32_t d = (int32_t)j - (int32_t)i;
    if (d < -(int32_t)w || d > (int32_t)w) return UINT16_MAX;
    return outIC[(size_t)i*stride + (uint32_t)(d + (int32_t)w)];
}

static inline uint16_t band2_get_notchecked(const uint16_t* outIC, uint32_t i, uint32_t j) {
    int32_t d = (int32_t)j - (int32_t)i;
    return outIC[(size_t)i*stride + (uint32_t)(d + (int32_t)w)];
}

// outIC layout: outIC[i*stride + (delta+w)] where delta=j-i in [-w..w]
#if 0
static inline void band2_distances_avx2_u16_xyz_v4(
	const uint16_t* __restrict xyz,   // [x...][y...][z...], each length N
	uint32_t N,
	uint16_t* __restrict outIC         // length N*(2*(M-1)+1)
){
	assert(xyz && outIC);
	if (N < 2 || M <= 1) return;

	const uint32_t w = M - 1;
	const uint32_t stride = 2*w + 1;

	const uint16_t* __restrict x = xyz;
	const uint16_t* __restrict y = xyz + N;
	const uint16_t* __restrict z = xyz + 2u*N;

	const __m256i i10000 = _mm256_set1_epi32(10000);
	const __m256i izero  = _mm256_setzero_si256();
	const __m256i i65535 = _mm256_set1_epi32(65535);

	for (uint32_t i = 0; i < N; ++i) {
		uint16_t* row = outIC + (size_t)i * stride;

		// Defensive: fill diagonal deterministically (distance 0 => IC=10000)
		row[w] = 10000;

		const __m256i xi = _mm256_set1_epi32((int32_t)x[i]);
		const __m256i yi = _mm256_set1_epi32((int32_t)y[i]);
		const __m256i zi = _mm256_set1_epi32((int32_t)z[i]);

		// Compute j range: [i-w, i+w] intersect [0,N-1]
		uint32_t j0 = (i > w) ? (i - w) : 0;
		uint32_t j1 = std::min(N - 1, i + w);

		// Compute forward j increasing; store into row at offset (j-i+w)
		uint32_t j = j0;

		// SIMD in chunks of 8 j's
		for (; j + 7 <= j1; j += 8) {
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

			// deterministic round-to-nearest
			fs = _mm256_round_ps(fs, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
			__m256i s_int = _mm256_cvttps_epi32(fs);

			__m256i dic = _mm256_add_epi32(s_int, i10000);
			dic = _mm256_max_epi32(dic, izero);
			dic = _mm256_min_epi32(dic, i65535);

			// pack 8x i32 -> 8x u16 contiguous
			__m128i dic_lo = _mm256_castsi256_si128(dic);
			__m128i dic_hi = _mm256_extracti128_si256(dic, 1);
			__m128i out8u16 = _mm_packus_epi32(dic_lo, dic_hi);

			// store at row[(j-i+w) .. +7]
			uint32_t off = (uint32_t)((int32_t)j - (int32_t)i + (int32_t)w);
			_mm_storeu_si128((__m128i*)(row + off), out8u16);
		}

		// scalar tail
		for (; j <= j1; ++j) {
			int32_t dx = (int32_t)x[j] - (int32_t)x[i];
			int32_t dy = (int32_t)y[j] - (int32_t)y[i];
			int32_t dz = (int32_t)z[j] - (int32_t)z[i];
			int32_t d2 = dx*dx + dy*dy + dz*dz;

			float s = std::sqrt((float)d2);
			uint32_t si = (uint32_t)std::floor(s + 0.5f);
			uint32_t dic = si + 10000u;
			if (dic > 65535u) dic = 65535u;

			uint32_t off = (uint32_t)((int32_t)j - (int32_t)i + (int32_t)w);
			row[off] = (uint16_t)dic;
		}
	}
}
#endif 

extern const uint32_t M; // you have: const uint32_t M = 100;

static inline void band2_distances_avx2_u16_xyz_v5(
    const uint16_t* __restrict xyz,     // length 3*N: [x...][y...][z...]
    uint32_t N,
    uint16_t* __restrict band,          // length N*stride, stride=2*(M-1)+1
    uint16_t* __restrict prev_neighbor, // length N, UINT16_MAX if none
    uint16_t* __restrict next_neighbor  // length N, UINT16_MAX if none
) {
    //const uint32_t w = M - 1;
    //const uint32_t stride = 2 * w + 1;

    const uint16_t* __restrict x = xyz;
    const uint16_t* __restrict y = xyz + N;
    const uint16_t* __restrict z = xyz + 2u * N;

    const __m256i i10000 = _mm256_set1_epi32(10000);
    const __m256i izero  = _mm256_setzero_si256();
    const __m256i i65535 = _mm256_set1_epi32(65535);

    alignas(32) uint16_t tmp_dic[8]; // for cheap horizontal scan

    for (uint32_t i = 0; i < N; ++i) {
        uint16_t* __restrict row = band + (size_t)i * stride;

        // band window for this row
        const uint32_t j0 = (i > w) ? (i - w) : 0;
        const uint32_t j1 = std::min(N - 1, i + w);

        // Fill diagonal deterministically (distance 0 => IC 10000)
        row[w] = 10000;

        // Initialize neighbors
        uint16_t best_prev_d = UINT16_MAX;
        uint16_t best_next_d = UINT16_MAX;
        uint16_t best_prev_j = UINT16_MAX;
        uint16_t best_next_j = UINT16_MAX;

        // Broadcast xi,yi,zi
        const __m256i xi = _mm256_set1_epi32((int32_t)x[i]);
        const __m256i yi = _mm256_set1_epi32((int32_t)y[i]);
        const __m256i zi = _mm256_set1_epi32((int32_t)z[i]);

        uint32_t j = j0;

        // SIMD over 8 j's
        for (; j + 7 <= j1; j += 8) {
            // load 8 coords
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

            __m256i dic32 = _mm256_add_epi32(s_int, i10000);
            dic32 = _mm256_max_epi32(dic32, izero);
            dic32 = _mm256_min_epi32(dic32, i65535);

            // pack 8x i32 -> 8x u16
            __m128i lo = _mm256_castsi256_si128(dic32);
            __m128i hi = _mm256_extracti128_si256(dic32, 1);
            __m128i dic16 = _mm_packus_epi32(lo, hi);

            // store distances into band row
            uint32_t off = (uint32_t)((int32_t)j - (int32_t)i + (int32_t)w);
            _mm_storeu_si128((__m128i*)(row + off), dic16);

            // update prev/next neighbors
            _mm_store_si128((__m128i*)tmp_dic, dic16);

            // lanes correspond to jj = j + lane
            // skip diagonal lane where jj==i
            for (uint32_t lane = 0; lane < 8; ++lane) {
                uint32_t jj = j + lane;
                if (jj == i) continue;
                uint16_t d = tmp_dic[lane];

                if (jj < i) {
                    if (d < best_prev_d) {
                        best_prev_d = d;
                        best_prev_j = (uint16_t)jj;
                    }
                } else { // jj > i
                    if (d < best_next_d) {
                        best_next_d = d;
                        best_next_j = (uint16_t)jj;
                    }
                }
            }
        }

        // scalar tail
        for (; j <= j1; ++j) {
            int32_t dx = (int32_t)x[j] - (int32_t)x[i];
            int32_t dy = (int32_t)y[j] - (int32_t)y[i];
            int32_t dz = (int32_t)z[j] - (int32_t)z[i];
            int32_t d2 = dx*dx + dy*dy + dz*dz;

            float s = std::sqrt((float)d2);
            uint32_t si = (uint32_t)std::floor(s + 0.5f);
            uint32_t dic = si + 10000u;
            if (dic > 65535u) dic = 65535u;

            uint32_t off = (uint32_t)((int32_t)j - (int32_t)i + (int32_t)w);
            row[off] = (uint16_t)dic;

            if (j == i) continue;
            if (j < i) {
                if ((uint16_t)dic < best_prev_d) {
                    best_prev_d = (uint16_t)dic;
                    best_prev_j = (uint16_t)j;
                }
            } else {
                if ((uint16_t)dic < best_next_d) {
                    best_next_d = (uint16_t)dic;
                    best_next_j = (uint16_t)j;
                }
            }
        }

        prev_neighbor[i] = best_prev_j;
        next_neighbor[i] = best_next_j;
    }
}