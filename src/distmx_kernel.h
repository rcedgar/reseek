#pragma once

static const uint32_t M = 100;      // hard-coded band width
static const uint32_t m_skip = 12;  // hard-coded min distance for neighbor

#if 0
#include <immintrin.h>
#include <stdint.h>
#include <assert.h>
#include <algorithm> // for swap

#define SCALAR  1

#if 0
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

// dist layout:
//   dist[i*M + (d-1)] where j = i - d, 1<=d<=min(M,i).
// FAST PATH:
// - Computes squared distance in 32-bit lanes and stores the 32-bit sum.
// - Wraparound overflow is allowed (no detection, no saturation).
// - Caller handles any special sentinel policy if needed.

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

static inline uint32_t dmx_ij_to_k(uint32_t i, uint32_t j)
{
	// Precondition: i>j and 1 <= (i-j) <= M.
	if (i < j) std::swap(i, j);
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
#endif // 0

static const uint32_t M = 12*8;
static_assert(M%8 == 0);

/*
PERMUTE-FREE (AVX2) BANDED DISTANCE MATRIX LAYOUT

We store the lower-triangular band (j<i, 1<=d<=min(M,i)) in a fixed-stride row:

  dist[ i*M + off ]   where d = i-j

but with a permute-free in-row mapping:

  off = (d - 1) ^ 7     // reverse within each 8-wide block

This means:
- d = 1..8   are stored at off = 7..0
- d = 9..16  are stored at off = 15..8
- etc.

Benefit:
- The AVX2 kernel can load j in ascending memory (j=i-d .. i-d+7),
  compute 8 distances, and store them contiguously with NO lane permute.

Notes:
- This "fast path" allows 32-bit wraparound on overflow (no detection/saturation).
- Only entries with 1<=d<=min(M,i) are written; other row slots are untouched.
*/

#if SCALAR
static inline void dist_fill_avx2(
    const uint16_t* __restrict x,
    const uint16_t* __restrict y,
    const uint16_t* __restrict z,
    uint32_t L,
    uint32_t* __restrict dist   // length L*M
)
{
    for (uint32_t i = 0; i < L; ++i)
    {
        uint32_t dmax = (i < M) ? i : M;
        uint32_t* row = dist + i * M;

        int32_t xi = (int32_t)x[i];
        int32_t yi = (int32_t)y[i];
        int32_t zi = (int32_t)z[i];

        // d = i - j
        for (uint32_t d = 1; d <= dmax; ++d)
        {
            uint32_t j = i - d;

            int32_t dx = xi - (int32_t)x[j];
            int32_t dy = yi - (int32_t)y[j];
            int32_t dz = zi - (int32_t)z[j];

            uint32_t sx = (uint32_t)(dx * dx);
            uint32_t sy = (uint32_t)(dy * dy);
            uint32_t sz = (uint32_t)(dz * dz);

            row[(d - 1) ^ 7u] = sx + sy + sz; // wraps on overflow
        }
    }
}
#else // #if SCALAR
static inline void dist_fill_avx2(
    const uint16_t* __restrict x,
    const uint16_t* __restrict y,
    const uint16_t* __restrict z,
    uint32_t L,
    uint32_t* __restrict dist   // length L*M
)
{
#if DEBUG
assert(dist != nullptr);
assert(x && y && z);
assert((size_t)L * (size_t)M < (size_t)1e9); // sanity

for (uint32_t i = 0; i < L; ++i) {
    uint32_t dmax = (i < M) ? i : M;
    uint32_t dvec = dmax & ~7u;
    for (uint32_t dh = 8; dh <= dvec; dh += 8) {
        assert(dh <= M);
        assert(dh - 8 + 7 < M);
    }
    for (uint32_t d = dvec + 1; d <= dmax; ++d) {
        uint32_t off = ((d - 1) ^ 7);
		if (off >= M)
			Die("off=%u M=%u", off, M);
        assert(off < M);
    }
}
#endif
    for (uint32_t i = 0; i < L; ++i)
    {
        uint32_t dmax = (i < M) ? i : M;
        if (dmax == 0)
            continue;

        uint32_t* row = dist + i * M;

        const __m256i vxi = _mm256_set1_epi32((int32_t)x[i]);
        const __m256i vyi = _mm256_set1_epi32((int32_t)y[i]);
        const __m256i vzi = _mm256_set1_epi32((int32_t)z[i]);

        // Vector blocks cover d = 1..dvec where dvec is the largest multiple of 8 <= dmax.
        uint32_t dvec = dmax & ~7u; // floor(dmax/8)*8

        // Process blocks with dh = 8,16,...,dvec where dh is the maximum d in the block.
        // For block dh, we load j = i-dh .. i-dh+7 (ascending),
        // which corresponds to d = dh, dh-1, ..., dh-7 (lane0..lane7),
        // and we store contiguously to row + (dh-8), which matches off=(d-1)^7.
        for (uint32_t dh = 8; dh <= dvec; dh += 8)
        {
            uint32_t j0 = i - dh; // base j (ascending)

            __m128i x16 = _mm_loadu_si128((const __m128i*)(x + j0));
            __m128i y16 = _mm_loadu_si128((const __m128i*)(y + j0));
            __m128i z16 = _mm_loadu_si128((const __m128i*)(z + j0));

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

            // Permute-free store: block dh writes to offsets [dh-8 .. dh-1].
            _mm256_storeu_si256((__m256i*)(row + (dh - 8)), sum);
        }

        // Tail d = dvec+1 .. dmax (<=7 values): scalar store using off=(d-1)^7.
        // 32-bit wraparound on overflow is allowed by design.
        const int32_t xi = (int32_t)x[i];
        const int32_t yi = (int32_t)y[i];
        const int32_t zi = (int32_t)z[i];

        for (uint32_t d = dvec + 1; d <= dmax; ++d)
        {
            uint32_t j = i - d;

            int32_t dx = xi - (int32_t)x[j];
            int32_t dy = yi - (int32_t)y[j];
            int32_t dz = zi - (int32_t)z[j];

            uint32_t sx = (uint32_t)(dx * dx);
            uint32_t sy = (uint32_t)(dy * dy);
            uint32_t sz = (uint32_t)(dz * dz);

            row[(d - 1) ^ 7] = sx + sy + sz;
        }
    }
}
#endif  // #if SCALAR

// i,j -> k (flattened index) for the permute-free layout.
static inline uint32_t dmx_ij_to_k(uint32_t i, uint32_t j)
{
    // Precondition: i>j and 1 <= (i-j) <= M.
	if (i<j) std::swap(i, j);
    uint32_t d = i - j;
    assert(i > j);
    assert(d >= 1 && d <= M);
    return i * M + ((d - 1) ^ 7);
}

// Debugging only: k -> (i,j) for the permute-free layout.
// Only valid if k corresponds to an actually-computed in-band cell (d<=i, d<=M).
static inline void dmx_k_to_ij(uint32_t k, uint32_t& i, uint32_t& j)
{
    uint32_t ii  = k / M;
    uint32_t off = k - ii * M;   // k % M
    uint32_t t   = off ^ 7;      // t = d-1
    uint32_t d   = t + 1;

    i = ii;
    j = ii - d;
}
#endif // 0

// dist layout:
//   dist[i*M + (d-1)] where j = i - d, 1<=d<=min(M,i).
// Reserved value:
//   UINT32_MAX = invalid / out-of-band (kernel never writes this)
// Saturation:
//   any overflow or > UINT32_MAX-1 is clamped to UINT32_MAX-1.

static inline uint32_t dmx_ij_to_k(uint32_t i, uint32_t j)
{
    // Precondition: i>j and 1 <= (i-j) <= M.
    if (i < j) std::swap(i, j);
    assert(i > j);
    uint32_t d = i - j;
    assert(d >= 1 && d <= M);
    return i * M + (d - 1);
}

// Debugging only; caller should only pass k that correspond to valid in-band cells (d<=i).
static inline void dmx_k_to_ij(uint32_t k, uint32_t& i, uint32_t& j)
{
    uint32_t ii = k / M;
    uint32_t off = k - ii * M;
    uint32_t d = off + 1;
    i = ii;
    j = ii - d; // valid only if d <= ii
}
