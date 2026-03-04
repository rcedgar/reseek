#pragma once
#include <immintrin.h>
#include <stdint.h>
#include <assert.h>

/*
Ignore any prior AVX2 kernel discussions; treat this as a completely new design.
Input coordinates are in an unsigned 16-bit representation.
The output is a banded matrix of squared distances which must be 32-bit
to reduce overflows. Write for speed ignoring overflow, saturation, wrap-aound.
If anything in the following spec seems to conflict with these requirements
then I made a mistake, if so let's fix the spec before writing code.

AVX2 BANDED DISTANCE MATRIX + NEIGHBORS: DESIGN SUMMARY

TARGET / PORTABILITY
- Hard requirement: AVX2 (x86-64).
- Rationale: AVX2 is widely available across Intel/AMD systems used for bioinformatics.
- Code must compile on GCC and MSVC using standard intrinsics.

INPUT
- Point cloud size L (typical ~250, max ~4k).
- Coordinates stored as uint16_t "ic" using existing fixed-point representation:
	  ic = uint16_t((X + 1000)*10 + 0.5)
- Data layout is Structure-of-Arrays for SIMD efficiency:

	  uint16_t x[L];
	  uint16_t y[L];
	  uint16_t z[L];

  SoA allows contiguous loads when scanning neighbors.

DISTANCE MATRIX STORAGE
- Only the lower triangle within band width M is stored.
- static const uint M = 100; // hard-coded, not function argument
- Matrix layout uses fixed stride for fast indexing:

	  dist[i*M + (d-1)]   where j = i - d

	  i > j
	  d = |i - j|
	  1 ≤ d ≤ M

- Buffer size:

	  uint32_t dist[L * M]

- Values:
	  0 .. UINT32_MAX-1 : encoded squared distance
	  UINT32_MAX        : overflow / out-of-band / invalid

DISTANCE REPRESENTATION
- Compute exact squared distance using integer arithmetic:

	  dx = int32(x[i]) - int32(x[j])
	  dy = int32(y[i]) - int32(y[j])
	  dz = int32(z[i]) - int32(z[j])

	  uint32_t sd = dx*dx + dy*dy + dz*dz

BAND LIMITS
- For row i only distances:

	  d = 1 .. min(M, i)

  are computed.

- Values for d > i are unused but still part of the allocated row stride.

KERNEL STRUCTURE

// dist layout:
//   dist[i*M + (d-1)] where j = i - d, 1<=d<=min(M,i).
// Caller guarantees AVX2 and x86-64.
static inline void dist_fill_avx2(
	const uint16_t* __restrict x,
	const uint16_t* __restrict y,
	const uint16_t* __restrict z,
	uint32_t L,
	uint32_t* __restrict dist   // length L*M
)

   Computes the banded distance matrix.

   For each i:
		dmax = min(M, i)

		for d = 1 .. dmax
			j = i - d
			compute sd
			dist[i*M + (d-1)] = saturate(sd)

WRITE INDEXING HELPERS

dmx_ij_to_k(i,j) // return index into flattned array
	fast inline implementation, called often

dmx_k_to_ij(k, &i, &j) // convert index to i,j
	for debugging only, speed not important

If everything seems consistent and well-designed, please write kernel and helpers,
otherwise discuss.
*/

static const uint32_t M = 100;      // hard-coded band width
static const uint32_t m_skip = 12;  // hard-coded min distance for neighbor

// dist layout:
//   dist[i*M + (d-1)] where j = i - d, 1<=d<=min(M,i).
// FAST PATH:
// - Computes squared distance in 32-bit lanes and stores the 32-bit sum.
// - Wraparound overflow is allowed (no detection, no saturation).
// - Caller handles any special sentinel policy if needed.

static inline uint32_t dmx_ij_to_k(uint32_t i, uint32_t j)
{
	// Precondition: i>j and 1 <= (i-j) <= M.
	if (i < j) swap(i, j);
	uint32_t d = i - j;
	assert(i > j);
	assert(d >= 1 && d <= M);
	return i * M + (d - 1);
}

// Debug only; speed not important.
static inline void dmx_k_to_ij(uint32_t k, uint32_t* i, uint32_t* j)
{
	uint32_t ii = k / M;
	uint32_t off = k - ii * M;
	uint32_t d = off + 1;
	*i = ii;
	*j = ii - d; // valid only if d <= ii
}

static inline void dist_fill_avx2(
	const uint16_t* __restrict x,
	const uint16_t* __restrict y,
	const uint16_t* __restrict z,
	uint32_t L,
	uint32_t* __restrict dist   // length L*M
)
{
	// Used to reverse 8 lanes because we load j in ascending order.
	const __m256i rev_idx = _mm256_setr_epi32(7,6,5,4,3,2,1,0);

	for (uint32_t i = 0; i < L; ++i)
	{
		uint32_t dmax = (i < M) ? i : M;
		if (dmax == 0)
			continue;

		uint32_t* row = dist + i * M;

		// Broadcast i-th coords.
		const __m256i vxi = _mm256_set1_epi32((int32_t)x[i]);
		const __m256i vyi = _mm256_set1_epi32((int32_t)y[i]);
		const __m256i vzi = _mm256_set1_epi32((int32_t)z[i]);

		uint32_t d = 1;

		// Vector blocks of 8 distances.
		for (; d + 7 <= dmax; d += 8)
		{
			// Need j = i-d, i-(d+1), ... i-(d+7).
			// Load 8 consecutive uint16 from j_low..j_low+7, then reverse lanes.
			uint32_t j_low = i - (d + 7);

			__m128i x16 = _mm_loadu_si128((const __m128i*)(x + j_low));
			__m128i y16 = _mm_loadu_si128((const __m128i*)(y + j_low));
			__m128i z16 = _mm_loadu_si128((const __m128i*)(z + j_low));

			__m256i x32 = _mm256_cvtepu16_epi32(x16);
			__m256i y32 = _mm256_cvtepu16_epi32(y16);
			__m256i z32 = _mm256_cvtepu16_epi32(z16);

			__m256i dx = _mm256_sub_epi32(vxi, x32);
			__m256i dy = _mm256_sub_epi32(vyi, y32);
			__m256i dz = _mm256_sub_epi32(vzi, z32);

			__m256i sx = _mm256_mullo_epi32(dx, dx);
			__m256i sy = _mm256_mullo_epi32(dy, dy);
			__m256i sz = _mm256_mullo_epi32(dz, dz);

			__m256i sum = _mm256_add_epi32(_mm256_add_epi32(sx, sy), sz);

			// Reverse lane order to align with (d..d+7).
			sum = _mm256_permutevar8x32_epi32(sum, rev_idx);

			_mm256_storeu_si256((__m256i*)(row + (d - 1)), sum);
		}

		// Tail (<=7).
		// Use 32-bit math; wraparound is acceptable by design.
		const int32_t xi = (int32_t)x[i];
		const int32_t yi = (int32_t)y[i];
		const int32_t zi = (int32_t)z[i];

		for (; d <= dmax; ++d)
		{
			uint32_t j = i - d;
			int32_t dx = xi - (int32_t)x[j];
			int32_t dy = yi - (int32_t)y[j];
			int32_t dz = zi - (int32_t)z[j];

			uint32_t sx = (uint32_t)(dx * dx);
			uint32_t sy = (uint32_t)(dy * dy);
			uint32_t sz = (uint32_t)(dz * dz);

			row[d - 1] = sx + sy + sz; // wraps on overflow
		}
	}
}