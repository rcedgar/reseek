#pragma once

#include <immintrin.h>
#include <stdint.h>
#include <string.h>   // memset

#ifndef asserta
  #include <assert.h>
  #define asserta(x) assert(x)
#endif

/*
AVX2 BANDED DISTANCE MATRIX + NEIGHBORS: DESIGN SUMMARY

TARGET / PORTABILITY
- Hard requirement: AVX2 (x86-64).
- Rationale: AVX2 is widely available across Intel/AMD systems used for bioinformatics.
  AVX-512 is inconsistent across CPUs and showed no advantage on Zen4 (Threadripper 7980X).
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
- Matrix layout uses fixed stride for fast indexing:

	  dist[i*M + (d-1)]   where j = i - d

	  i > j
	  d = |i - j|
	  1 ≤ d ≤ M

- Buffer size:

	  uint16_t dist[L * M]

- Values:
	  0 .. UINT16_MAX-1 : encoded squared distance
	  UINT16_MAX        : out-of-band / invalid

DISTANCE REPRESENTATION
- Compute exact squared distance using integer arithmetic:

	  dx = int32(x[i]) - int32(x[j])
	  dy = int32(y[i]) - int32(y[j])
	  dz = int32(z[i]) - int32(z[j])

	  sd = dx*dx + dy*dy + dz*dz      // uint32 intermediate

- Store result in uint16_t with saturation:

	  dist_code = min(sd, UINT16_MAX-1)

- UINT16_MAX is reserved exclusively for "out of band".

Rationale
- Using uint16 storage reduces memory bandwidth and improves cache locality.
- Exact squared distances are preserved for the most important small values.
- Very large distances saturate but those are not relevant for neighbor selection.
- Integer comparison preserves ordering except when saturation occurs.

BAND LIMITS
- For row i only distances:

	  d = 1 .. min(M, i)

  are computed.

- Values for d > i are unused but still part of the allocated row stride.

NEIGHBOR DEFINITIONS
Two nearest-neighbor vectors are required:

1) predecessor neighbor

	  prev[i] =
		  argmin_{j < i, (i-j) > m_skip, (i-j) ≤ M} dist(i,j)

2) successor neighbor

	  next[i] =
		  argmin_{j > i, (j-i) > m_skip, (j-i) ≤ M} dist(i,j)

- If no valid neighbor exists, store UINT16_MAX.

Typical parameters:
	  M ≈ 100
	  m_skip ≈ 12

KERNEL STRUCTURE

1) dist_fill_avx2()

   Computes the banded distance matrix.

   For each i:
		dmax = min(M, i)

		for d = 1 .. dmax
			j = i - d
			compute sd
			dist[i*M + (d-1)] = saturate(sd)

   SIMD strategy:
   - Process 16 neighbor points per iteration using AVX2.
   - Broadcast (x[i],y[i],z[i]).
   - Load contiguous x[j..], y[j..], z[j..].
   - Compute dx,dy,dz in epi16 then widen to epi32.
   - Multiply/add to obtain squared distances.
   - Saturate to uint16 and store.

   This kernel is the main performance hotspot and benefits strongly from SIMD.

2) prev_from_dist()

   Computes predecessor neighbors.

   For each i:
	   scan row dist[i*M + d-1] for d = m_skip+1 .. min(M,i)
	   choose minimum distance

   This is a contiguous row reduction and can optionally be vectorized.

3) next_from_dist()

   Computes successor neighbors.

   Iterate through matrix entries:

	   for i = 0..L-1
		   for d = m_skip+1 .. min(M,i)
			   j = i - d
			   code = dist[i*M + (d-1)]
			   update best successor for j

   This pass involves scattered writes and is therefore not SIMD-friendly.

FUSION DECISION

Distance computation is **not fused** with neighbor updates.

Reason:
- Distance computation is highly SIMD-friendly.
- Neighbor updates require per-lane scalar work and scattered memory updates.
- Fusing them destroys the SIMD speedup (empirically confirmed).

Separating kernels allows the distance kernel to run ~3× faster on AVX2.

INDEXING HELPERS

(i,j) → k

	  hi = max(i,j)
	  lo = min(i,j)
	  d  = hi - lo

	  if d == 0 or d > M → invalid

	  k = hi*M + (d-1)

k → (i,j) (debug only)

	  i = k / M
	  d = (k % M) + 1
	  j = i - d

NEXT IMPLEMENTATION STEP

Implement three kernels:

	  dist_fill_avx2()
	  prev_from_dist()
	  next_from_dist()
*/

// Reserved codes
static constexpr uint16_t DIST_OOB   = UINT16_MAX;     // out-of-band / invalid
static constexpr uint16_t DIST_SAT   = UINT16_MAX - 1; // saturation cap (largest valid code)

// Reverse 16x uint16 lanes in a __m256i (full 256-bit reversal).
static inline __m256i reverse_u16x16(__m256i v)
{
	// swap 128-bit halves
	v = _mm256_permute2x128_si256(v, v, 0x01);

	// reverse bytes within each 128-bit half in 16-bit chunks:
	// [0..15] u16 -> [15..0]
	const __m256i shuf = _mm256_setr_epi8(
		14,15, 12,13, 10,11,  8, 9,  6, 7,  4, 5,  2, 3,  0, 1,
		14,15, 12,13, 10,11,  8, 9,  6, 7,  4, 5,  2, 3,  0, 1
	);
	return _mm256_shuffle_epi8(v, shuf);
}

// dist layout:
//   dist[i*M + (d-1)] where j = i - d, 1<=d<=min(M,i).
// Caller guarantees AVX2 and x86-64.
static inline void dist_fill_avx2(
	const uint16_t* __restrict x,
	const uint16_t* __restrict y,
	const uint16_t* __restrict z,
	uint32_t L,
	uint32_t M,
	uint16_t* __restrict dist   // length L*M
)
{
#if DEBUG
	// Optional: pre-fill with OOB so unused (d>i) entries are consistently invalid.
	// This can help debug and simplifies downstream if code ever touches unused cells.
	for (uint32_t i = 0; i < L; ++i) {
		uint16_t* row = dist + (size_t)i * M;
		for (uint32_t k = 0; k < M; ++k) row[k] = DIST_OOB;
	}
#endif
	const __m256i vsat = _mm256_set1_epi32((int)DIST_SAT);

	for (uint32_t i = 0; i < L; ++i)
	{
		const uint32_t dmax = (i < M) ? i : M;
		if (dmax == 0) continue;

		const __m256i xi16 = _mm256_set1_epi16((short)x[i]);
		const __m256i yi16 = _mm256_set1_epi16((short)y[i]);
		const __m256i zi16 = _mm256_set1_epi16((short)z[i]);

		uint16_t* row = dist + (size_t)i * M;

		// SIMD blocks: we load contiguous j in ascending order, but d = i - j
		// decreases across lanes. Compute, then reverse the packed u16 results
		// so we can store to row[(d-1) .. (d-16)] correctly as contiguous.
		uint32_t d = 1;
		for (; d + 15 <= dmax; d += 16)
		{
			// We want d..d+15. Corresponding j are i-d .. i-(d+15) (descending).
			// Load contiguous ascending block [j_lo .. j_hi] where:
	
			const uint32_t j_lo = i - (d + 15);   // smallest index
			// const uint32_t j_hi = i - d;       // largest index

			// Load 16x u16 for each coordinate
			__m256i xj16 = _mm256_loadu_si256((const __m256i*)(x + j_lo));
			__m256i yj16 = _mm256_loadu_si256((const __m256i*)(y + j_lo));
			__m256i zj16 = _mm256_loadu_si256((const __m256i*)(z + j_lo));

			// dx/dy/dz as signed 16-bit (assumes |diff| fits in int16; true for your ic scaling)
			__m256i dx16 = _mm256_sub_epi16(xi16, xj16);
			__m256i dy16 = _mm256_sub_epi16(yi16, yj16);
			__m256i dz16 = _mm256_sub_epi16(zi16, zj16);

			// Widen low/high 8 lanes to 32-bit
			__m256i dx_lo = _mm256_cvtepi16_epi32(_mm256_castsi256_si128(dx16));
			__m256i dx_hi = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(dx16, 1));
			__m256i dy_lo = _mm256_cvtepi16_epi32(_mm256_castsi256_si128(dy16));
			__m256i dy_hi = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(dy16, 1));
			__m256i dz_lo = _mm256_cvtepi16_epi32(_mm256_castsi256_si128(dz16));
			__m256i dz_hi = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(dz16, 1));

			// sd = dx*dx + dy*dy + dz*dz (fits in uint32)
			__m256i sd_lo = _mm256_mullo_epi32(dx_lo, dx_lo);
			__m256i sd_hi = _mm256_mullo_epi32(dx_hi, dx_hi);
			sd_lo = _mm256_add_epi32(sd_lo, _mm256_mullo_epi32(dy_lo, dy_lo));
			sd_hi = _mm256_add_epi32(sd_hi, _mm256_mullo_epi32(dy_hi, dy_hi));
			sd_lo = _mm256_add_epi32(sd_lo, _mm256_mullo_epi32(dz_lo, dz_lo));
			sd_hi = _mm256_add_epi32(sd_hi, _mm256_mullo_epi32(dz_hi, dz_hi));

			// clamp to DIST_SAT (65534) so we never emit UINT16_MAX (reserved OOB)
			sd_lo = _mm256_min_epu32(sd_lo, vsat);
			sd_hi = _mm256_min_epu32(sd_hi, vsat);

			// Pack to 16x u16 (unsigned-saturating); values are <= 65534 so safe.
			__m256i packed = _mm256_packus_epi32(sd_lo, sd_hi);
			// packus produces lane order grouped by 128-bit halves; fix ordering:
			packed = _mm256_permute4x64_epi64(packed, 0xD8);

			// Our loaded j block corresponds to d+15 .. d (descending),
			// so reverse packed u16 to get d .. d+15.
			packed = reverse_u16x16(packed);

			_mm256_storeu_si256((__m256i*)(row + (d - 1)), packed);
		}

		// Tail (scalar)
		for (; d <= dmax; ++d)
		{
			const uint32_t j = i - d;
			const int32_t dx = (int32_t)x[i] - (int32_t)x[j];
			const int32_t dy = (int32_t)y[i] - (int32_t)y[j];
			const int32_t dz = (int32_t)z[i] - (int32_t)z[j];
			const uint32_t sd = (uint32_t)(dx*dx + dy*dy + dz*dz);
			row[d - 1] = (sd > (uint32_t)DIST_SAT) ? DIST_SAT : (uint16_t)sd;
		}
	}
}

// Convert matrix coordinates (i,j) → linear index k
// Layout:
//      dist[i*M + (d-1)]  where d = |i-j|, i>j
//
// checked: asserts if invalid or out of band.
static inline uint32_t dmx_ij_to_k_checked(
    uint32_t i,
    uint32_t j,
    uint32_t M)
{
    asserta(i != j);

    uint32_t hi = (i > j) ? i : j;
    uint32_t lo = (i > j) ? j : i;

    uint32_t d = hi - lo;
    //if (d == 0 || d > M)
    //    return UINT32_MAX;
	asserta(d != 0);
	asserta(d <= M);

    return hi*M + (d - 1);
}

static inline uint32_t dmx_ij_to_k(
    uint32_t i,
    uint32_t j,
    uint32_t M)
{
    uint32_t hi = (i > j) ? i : j;
    uint32_t lo = (i > j) ? j : i;
    uint32_t d = hi - lo;
    return hi*M + (d - 1);
}


// Convert linear index k → (i,j)
// Mainly for debugging / verification.
static inline void dmx_k_to_ij(
    uint32_t k,
    uint32_t M,
    uint32_t &i,
    uint32_t &j)
{
    i = k / M;
    uint32_t d = (k % M) + 1;

    j = i - d;
}