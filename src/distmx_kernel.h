#pragma once

/*
DISTMX SIMD KERNEL DESIGN SUMMARY
Highlights
- SoA coords, banded lower-triangle, fused neighbors

GOAL / WORKLOAD
- Build a symmetric distance matrix for a point cloud of size L (typical ~250, max ~4k).
- Only store a band of width M in cloud order (typical M~100): distances for |i-j| in [1..M].
- Matrix cells may be randomly accessed ~5-10x per DP/other kernel, so (i,j)->k indexing must be very cheap.
- A common downstream operation is nearest predecessor and nearest successor (by distance), excluding close
  indices |i-j| <= m_skip (typical m_skip~12). If no neighbor in band, use UINT16_MAX.
- NOTE m_skip applies ONLY to neighbor calculation, all off-diagonal values within band M must be stored.

INPUT COORDINATE FORMAT
- Coordinates are stored as uint16_t "ic" with offset +1000 and scale *10:
    ic = uint16_t((X + 1000)*10 + 0.5)
  This format is NOT changed by this design.

DISTANCE REPRESENTATION (squared distances)
- We store SQUARED Euclidean distances computed from ic deltas:
    dx = int32_t(x[i]) - int32_t(x[j])  // units: 0.1
    sd = dx*dx + dy*dy + dz*dz          // units: 0.01 (i.e., (0.1)^2)
- Stored type is uint16_t with two reserved values:
    0xFFFF = out-of-band / not computed
    0xFFFE = saturation ceiling for large distances (kept < 0xFFFF to distinguish from OOB)

- Encoding prioritizes accuracy / ordering for SMALL distances (Euclidean < 10.0),
  while tolerating coarser quantization for larger distances to extend dynamic range.

- Piecewise-linear monotone encoding from exact sd (uint32_t) to uint16_t code:
    Let T0 = 16383.  // exact region limit in sd units (0.01); corresponds to Euclid ~ sqrt(T0)/10 ~= 12.8
    Let T1 = 81915.  // end of mild-quantization region, T1 = T0 + (16383<<2)

    if (sd <= T0)            code = (uint16_t)sd;                          // exact ordering
    else if (sd <= T1)       code = T0 + ((sd - T0) >> 2);                 // step=4 in sd units
    else                     code = 32768 + ((sd - T1) >> 6);              // step=64 in sd units
                              code clamped to 0xFFFE                       // avoid 0xFFFF

  Properties:
  - Strictly monotone in sd (ties only from quantization), so neighbor ordering is preserved exactly for sd<=T0.
  - T0 chosen so Euclidean < 10.0 (sd < 10000) is fully inside the exact region.
  - Large distances are representable without early saturation; only very large sd values clamp to 0xFFFE.

COORDINATE STORAGE LAYOUT (FOR SIMD)
- Use Structure-of-Arrays (SoA):
    uint16_t x[L], y[L], z[L] (64B aligned preferred).
  Rationale: contiguous loads for j blocks; avoids de-interleaving/gathers; best for AVX2/AVX-512.

DISTANCE MATRIX STORAGE LAYOUT (FOR FAST RANDOM ACCESS)
- Store only LOWER band (j<i) with FIXED STRIDE = M entries per row i:
    dist[i*M + (d-1)] corresponds to j = i - d, with d in [1..M].
  Access:
    if (i==j) -> treat as 0 (not stored)
    else let (hi,lo) = (max(i,j), min(i,j)), d = hi-lo:
         if (d<1 || d>M) => OOB
         else k = hi*M + (d-1)
  Rationale: index is multiply+add only; ideal for random accesses and later kernels.

BAND / SKIP POLICY
- Only compute/store entries for d in [m_skip+1 .. min(M, i)] for each row i.
- Entries for d <= m_skip or d > min(M,i) are set to OOB (0xFFFF) for debugging and comparisons.
  (Performance option: skip storing OOB for d<=m_skip if caller never reads them.)

FUSED NEIGHBOR COMPUTATION (SAME PASS AS MATRIX FILL)
- While computing each (i>j) pair:
  - It is a predecessor candidate for i (prev[i]).
  - It is a successor candidate for j (next[j]).
- Maintain best distance seen so far per row i for prev, and per point j for next.
- Because AVX2 has no scatter store, next[] updates are done in a small per-lane scalar loop after SIMD
  distance computation; this is still typically faster than a second full matrix pass.

NEIGHBOR OUTPUTS
- prev[i] = argmin_j<i, (i-j)>m_skip, (i-j)<=M of encoded distance; UINT16_MAX if none.
- next[j] = argmin_i>j, (i-j)>m_skip, (i-j)<=M of encoded distance; UINT16_MAX if none.
- Tie-breaking (recommended for determinism): if encoded distances equal, prefer smaller |i-j|, then smaller index.

KERNEL VARIANTS
- Reference / clarity kernel (scalar): simplest loops, exact sd, uses same encoding and neighbor logic.
  Serves as correctness oracle.
- AVX2 kernel: baseline fast path, portable across MSVC and GCC/Clang with AVX2.
- Newer-hardware kernel: AVX-512 variant for higher throughput (wider vectors + masking), same output format.

VALIDATION SUPPORT
- Provide a comparer that checks dist[] and prev/next[] from two kernels and reports any mismatches.
- Since the encoding is pure integer and deterministic, mismatches are true bugs (not rounding noise).
*/

#define __AVX512F__
#define __AVX512BW__
#define __AVX2__

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <algorithm>
#include <limits>
#include <vector>

#if defined(_MSC_VER)
  #include <intrin.h>
#endif

#if defined(__AVX2__) || (defined(_MSC_VER) && defined(__AVX2__))
  #include <immintrin.h>
#endif

#if defined(__AVX512F__) && defined(__AVX512BW__)
  #include <immintrin.h>
#endif

namespace distmx {

// ============================
// Compile-time configuration
// ============================

// User-requested: static const, not function args.
static const uint32_t M      = 100;  // band width (max |i-j|)
static const uint32_t m_skip = 12;   // neighbor selection excludes |i-j| <= m_skip (distances still stored)

// Reserved distance codes
static const uint16_t DIST_OOB = 0xFFFF; // out-of-band / not computed
static const uint16_t DIST_SAT = 0xFFFE; // saturation ceiling (<DIST_OOB)

// Distance encoding thresholds (in squared ic-units)
// ic units: 0.1; squared units: 0.01
// Euclidean < 10.0 => sd < 10000 which is within exact region (sd <= T0).
static const uint32_t T0 = 16383;                 // exact region end  (~12.8 units)
static const uint32_t T1 = T0 + (16383u << 2);    // end of mild-quantization region (step=4)

// ============================
// Output struct
// ============================

struct Out {
  // dist: length L*M. dist[i*M + (d-1)] corresponds to j=i-d, for d=1..M.
  // For d>min(M,i) (i.e., j<0): DIST_OOB.
  // NOTE: distances for d<=m_skip are STILL stored; m_skip affects only neighbor selection.
  uint16_t* dist;   // [L*M]
  uint16_t* prev;   // [L] predecessor index (j<i) or UINT16_MAX
  uint16_t* next;   // [L] successor index (j>i) or UINT16_MAX
};

// ============================
// Helpers
// ============================

// Convert (i,j) to flat offset k in dist[].
// Layout: dist[i*M + (d-1)] where d = i-j, 1 <= d <= M, and i>j.
// Only lower triangle stored; symmetry handled internally.
//
// Returns:
//   - k (0 <= k < L*M) if in band
//   - SIZE_MAX if (i==j) or |i-j| > M
//
static inline uint32_t ij_to_k(uint32_t i, uint32_t j)
{
  assert(i != j);
  // ensure hi > lo
  uint32_t hi = (i > j) ? i : j;
  uint32_t lo = (i > j) ? j : i;
  uint32_t d = hi - lo;          // |i-j|
  assert(d != 0);
  assert(d < M);
  // row = hi, column offset = (d-1)
  return hi * M + (d - 1);
}

// Returns number of uint16_t elements required for dist[] buffer.
// Layout is fixed-stride lower band: L rows × M columns.
static inline uint32_t dist_buffer_length(uint32_t L)
{
  return size_t(L) * M;
}

static inline uint16_t encode_sd_to_u16(uint32_t sd) {
  // Piecewise-linear monotone encoding, prioritizing exactness for small distances.
  // Reserves 0xFFFF for OOB and clamps at 0xFFFE.
  uint32_t code;
  if (sd <= T0) {
    code = sd;
  } else if (sd <= T1) {
    code = T0 + ((sd - T0) >> 2);              // step=4
  } else {
    code = 32768u + ((sd - T1) >> 6);          // step=64
  }
  if (code >= (uint32_t)DIST_SAT) return DIST_SAT;
  return (uint16_t)code;
}

// Decode encoded uint16 squared-distance code to approximate Euclidean distance (float).
// Slow (uses sqrt); intended for debugging/logging.
// Units: same as original coordinates (before coord2ic scaling).
static inline float decode_code_to_distance(uint16_t code)
{
  if (code == DIST_OOB)
    return FLT_MAX;

  // Choose a representative sd (uint32) for this code.
  // For exact region, inversion is exact.
  // For quantized regions, choose midpoint of the bucket to reduce bias.
  uint32_t sd;

  if (code <= (uint16_t)T0) {
    sd = (uint32_t)code;
  } else if (code < 32768u) {
    // code = T0 + ((sd - T0) >> 2)
    // Let q = code - T0, so (sd - T0) in [q<<2, (q+1<<2)-1]
    uint32_t q = (uint32_t)code - T0;
    uint32_t lo = T0 + (q << 2);
    uint32_t hi = lo + 3;
    sd = (lo + hi) >> 1; // midpoint
  } else {
    // code = 32768 + ((sd - T1) >> 6)  (clamped at DIST_SAT)
    uint32_t q = (uint32_t)code - 32768u;
    uint32_t lo = T1 + (q << 6);
    uint32_t hi = lo + 63;
    sd = (lo + hi) >> 1; // midpoint
  }

  // Convert sd (ic^2 units; 0.01) -> Euclidean distance in coordinate units:
  // dx in ic is tenths, so sqrt(sd) gives tenths; multiply by 0.1.
  return std::sqrt((float)sd) * 0.1f;
}

static inline void init_outputs(uint32_t L, Out out) {
  std::fill(out.dist, out.dist + size_t(L) * M, DIST_OOB);
  std::fill(out.prev, out.prev + L, uint16_t(0xFFFF));
  std::fill(out.next, out.next + L, uint16_t(0xFFFF));
}

// Deterministic tie-break for neighbors:
// - lower encoded distance wins
// - if tie: smaller d (closer in cloud order) wins
// - if tie: smaller index wins
static inline bool better_candidate(uint16_t dist_new, uint32_t d_new, uint32_t idx_new,
                                    uint16_t dist_best, uint32_t d_best, uint32_t idx_best) {
  if (dist_new < dist_best) return true;
  if (dist_new > dist_best) return false;
  if (d_new < d_best) return true;
  if (d_new > d_best) return false;
  return idx_new < idx_best;
}

// Compare two outputs; returns true if identical. If not, writes first diffs to `buf` (optional).
// buf may be null; buf_len may be 0. max_reports limits output volume.
static inline bool compare_outputs(const Out& A, const Out& B,
                                   uint32_t L,
                                   char* buf, size_t buf_len,
                                   size_t max_reports = 50) {
  size_t reports = 0;
  size_t pos = 0;
  auto emit = [&](const char* msg) {
    if (!buf || buf_len == 0) return;
    size_t n = std::min(buf_len - pos, std::strlen(msg));
    if (n) { std::memcpy(buf + pos, msg, n); pos += n; }
    if (pos < buf_len) buf[pos] = '\0';
  };

  auto emit_line3 = [&](const char* tag, uint32_t i, uint32_t a, uint32_t b) {
    if (!buf || buf_len == 0) return;
    char line[160];
#if defined(_MSC_VER)
    _snprintf_s(line, sizeof(line), _TRUNCATE, "%s i=%u A=%u B=%u\n", tag, i, a, b);
#else
    std::snprintf(line, sizeof(line), "%s i=%u A=%u B=%u\n", tag, i, a, b);
#endif
    emit(line);
  };

  auto emit_line_k = [&](uint32_t k, uint32_t a, uint32_t b) {
    if (!buf || buf_len == 0) return;
    char line[160];
#if defined(_MSC_VER)
    _snprintf_s(line, sizeof(line), _TRUNCATE, "dist k=%u A=%u B=%u\n", k, a, b);
#else
    std::snprintf(line, sizeof(line), "dist k=%u A=%u B=%u\n", k, a, b);
#endif
    emit(line);
  };

  bool ok = true;

  for (uint32_t i = 0; i < L; ++i) {
    if (A.prev[i] != B.prev[i]) {
      ok = false;
      if (reports++ < max_reports) emit_line3("prev", i, A.prev[i], B.prev[i]);
      else break;
    }
    if (A.next[i] != B.next[i]) {
      ok = false;
      if (reports++ < max_reports) emit_line3("next", i, A.next[i], B.next[i]);
      else break;
    }
  }

  const uint32_t N = L * M;
  for (uint32_t k = 0; ok && k < N; ++k) {
    if (A.dist[k] != B.dist[k]) {
      ok = false;
      if (reports++ < max_reports) emit_line_k(k, A.dist[k], B.dist[k]);
      else break;
    }
  }

  return ok;
}

// ============================
// Variant 1: scalar reference
// ============================

static inline void kernel_scalar_ref(const uint16_t* x, const uint16_t* y, const uint16_t* z,
                                     uint32_t L, Out out) {
  init_outputs(L, out);

  // Successor tracking (encoded distance and d=i-j for tie-break).
  std::vector<uint16_t> best_next_d(L, DIST_OOB);
  std::vector<uint32_t> best_next_didx(L, std::numeric_limits<uint32_t>::max());

  for (uint32_t i = 0; i < L; ++i) {
    const uint32_t dmax = std::min<uint32_t>(M, i);

    uint16_t best_prev_dist = DIST_OOB;
    uint32_t best_prev_d = std::numeric_limits<uint32_t>::max();
    uint32_t best_prev_j = std::numeric_limits<uint32_t>::max();

    for (uint32_t d = 1; d <= dmax; ++d) {
      const uint32_t j = i - d;

      const int32_t dx = int32_t(x[i]) - int32_t(x[j]);
      const int32_t dy = int32_t(y[i]) - int32_t(y[j]);
      const int32_t dz = int32_t(z[i]) - int32_t(z[j]);

      const uint32_t sd = uint32_t(dx*dx + dy*dy + dz*dz);
      const uint16_t code = encode_sd_to_u16(sd);

      out.dist[size_t(i) * M + (d - 1)] = code;

      // Neighbor selection excludes close-by indices (d <= m_skip).
      if (d > m_skip) {
        // prev[i]
        if (out.prev[i] == 0xFFFF ||
            better_candidate(code, d, j, best_prev_dist, best_prev_d, best_prev_j)) {
          best_prev_dist = code;
          best_prev_d = d;
          best_prev_j = j;
          out.prev[i] = (uint16_t)j;
        }

        // next[j]
        if (out.next[j] == 0xFFFF ||
            better_candidate(code, d, i, best_next_d[j], best_next_didx[j], out.next[j])) {
          best_next_d[j] = code;
          best_next_didx[j] = d;
          out.next[j] = (uint16_t)i;
        }
      }
    }
  }
}

// ============================
// Variant 2: AVX2 baseline
// - Computes exact sd in SIMD (32-bit), stores encoded code per lane.
// - Neighbor selection is gated by (d > m_skip); distances always stored for d=1..dmax.
// ============================

#if defined(__AVX2__) || (defined(_MSC_VER) && defined(__AVX2__))

static inline void kernel_avx2(const uint16_t* x, const uint16_t* y, const uint16_t* z,
                               uint32_t L, Out out) {
  init_outputs(L, out);

  std::vector<uint16_t> best_next_d(L, DIST_OOB);
  std::vector<uint32_t> best_next_didx(L, std::numeric_limits<uint32_t>::max());

  alignas(32) uint32_t tmp_sd_lo[8];
  alignas(32) uint32_t tmp_sd_hi[8];

  for (uint32_t i = 0; i < L; ++i) {
    const uint32_t dmax = std::min<uint32_t>(M, i);

    uint16_t best_prev_dist = DIST_OOB;
    uint32_t best_prev_d = std::numeric_limits<uint32_t>::max();
    uint32_t best_prev_j = std::numeric_limits<uint32_t>::max();

    const __m256i vxi = _mm256_set1_epi16((int16_t)x[i]);
    const __m256i vyi = _mm256_set1_epi16((int16_t)y[i]);
    const __m256i vzi = _mm256_set1_epi16((int16_t)z[i]);

    uint32_t d = 1;

    // Process 16 distances per block (two halves of 8 widened to epi32).
    for (; d + 15 <= dmax; d += 16) {
      // j ranges [i-(d+15) .. i-d] increasing in memory
      const uint32_t j_low = i - (d + 15);

      __m256i vxj = _mm256_loadu_si256((const __m256i*)(x + j_low));
      __m256i vyj = _mm256_loadu_si256((const __m256i*)(y + j_low));
      __m256i vzj = _mm256_loadu_si256((const __m256i*)(z + j_low));

      __m256i vdx16 = _mm256_sub_epi16(vxi, vxj); // dx = xi - xj
      __m256i vdy16 = _mm256_sub_epi16(vyi, vyj);
      __m256i vdz16 = _mm256_sub_epi16(vzi, vzj);

      // low 8 lanes
      __m128i dx_lo = _mm256_castsi256_si128(vdx16);
      __m128i dy_lo = _mm256_castsi256_si128(vdy16);
      __m128i dz_lo = _mm256_castsi256_si128(vdz16);

      __m256i dx0 = _mm256_cvtepi16_epi32(dx_lo);
      __m256i dy0 = _mm256_cvtepi16_epi32(dy_lo);
      __m256i dz0 = _mm256_cvtepi16_epi32(dz_lo);

      __m256i s0 = _mm256_add_epi32(
                     _mm256_add_epi32(_mm256_mullo_epi32(dx0, dx0), _mm256_mullo_epi32(dy0, dy0)),
                     _mm256_mullo_epi32(dz0, dz0));
      _mm256_store_si256((__m256i*)tmp_sd_lo, s0);

      // high 8 lanes
      __m128i dx_hi = _mm256_extracti128_si256(vdx16, 1);
      __m128i dy_hi = _mm256_extracti128_si256(vdy16, 1);
      __m128i dz_hi = _mm256_extracti128_si256(vdz16, 1);

      __m256i dx1 = _mm256_cvtepi16_epi32(dx_hi);
      __m256i dy1 = _mm256_cvtepi16_epi32(dy_hi);
      __m256i dz1 = _mm256_cvtepi16_epi32(dz_hi);

      __m256i s1 = _mm256_add_epi32(
                     _mm256_add_epi32(_mm256_mullo_epi32(dx1, dx1), _mm256_mullo_epi32(dy1, dy1)),
                     _mm256_mullo_epi32(dz1, dz1));
      _mm256_store_si256((__m256i*)tmp_sd_hi, s1);

      // Store / neighbor updates in increasing d order.
      // Loaded j increasing corresponds to d decreasing; map carefully:
      // tmp_sd_lo lanes 0..7  => j=j_low..j_low+7   => d=(d+15)..(d+8)
      // tmp_sd_hi lanes 0..7  => j=j_low+8..j_low+15=> d=(d+7)..(d)
      // We want dt = d..d+15 (increasing), so reverse accordingly.
      for (uint32_t t = 0; t < 16; ++t) {
        const uint32_t dt = d + t;
        const uint32_t j  = i - dt;

        uint32_t sd = (t < 8) ? tmp_sd_hi[7 - t] : tmp_sd_lo[15 - t];
        const uint16_t code = encode_sd_to_u16(sd);

        out.dist[size_t(i) * M + (dt - 1)] = code;

        if (dt > m_skip) {
          // prev[i]
          if (out.prev[i] == 0xFFFF ||
              better_candidate(code, dt, j, best_prev_dist, best_prev_d, best_prev_j)) {
            best_prev_dist = code;
            best_prev_d = dt;
            best_prev_j = j;
            out.prev[i] = (uint16_t)j;
          }

          // next[j]
          if (out.next[j] == 0xFFFF ||
              better_candidate(code, dt, i, best_next_d[j], best_next_didx[j], out.next[j])) {
            best_next_d[j] = code;
            best_next_didx[j] = dt;
            out.next[j] = (uint16_t)i;
          }
        }
      }
    }

    // Tail scalar
    for (; d <= dmax; ++d) {
      const uint32_t j = i - d;
      const int32_t dx = int32_t(x[i]) - int32_t(x[j]);
      const int32_t dy = int32_t(y[i]) - int32_t(y[j]);
      const int32_t dz = int32_t(z[i]) - int32_t(z[j]);
      const uint32_t sd = uint32_t(dx*dx + dy*dy + dz*dz);
      const uint16_t code = encode_sd_to_u16(sd);

      out.dist[size_t(i) * M + (d - 1)] = code;

      if (d > m_skip) {
        if (out.prev[i] == 0xFFFF ||
            better_candidate(code, d, j, best_prev_dist, best_prev_d, best_prev_j)) {
          best_prev_dist = code;
          best_prev_d = d;
          best_prev_j = j;
          out.prev[i] = (uint16_t)j;
        }

        if (out.next[j] == 0xFFFF ||
            better_candidate(code, d, i, best_next_d[j], best_next_didx[j], out.next[j])) {
          best_next_d[j] = code;
          best_next_didx[j] = d;
          out.next[j] = (uint16_t)i;
        }
      }
    }
  }
}

#endif // AVX2

// ============================
// Variant 3: AVX-512BW (newer hardware)
// - Wider loads, same deterministic encoding.
// - Neighbor selection gated by (d > m_skip); distances always stored.
// ============================

#if defined(__AVX512F__) && defined(__AVX512BW__)

static inline void kernel_avx512bw(const uint16_t* x, const uint16_t* y, const uint16_t* z,
                                   uint32_t L, Out out) {
  init_outputs(L, out);

  std::vector<uint16_t> best_next_d(L, DIST_OOB);
  std::vector<uint32_t> best_next_didx(L, std::numeric_limits<uint32_t>::max());

  alignas(64) uint32_t tmp_sd_lo[16];
  alignas(64) uint32_t tmp_sd_hi[16];

  for (uint32_t i = 0; i < L; ++i) {
    const uint32_t dmax = std::min<uint32_t>(M, i);

    uint16_t best_prev_dist = DIST_OOB;
    uint32_t best_prev_d = std::numeric_limits<uint32_t>::max();
    uint32_t best_prev_j = std::numeric_limits<uint32_t>::max();

    const __m512i vxi = _mm512_set1_epi16((int16_t)x[i]);
    const __m512i vyi = _mm512_set1_epi16((int16_t)y[i]);
    const __m512i vzi = _mm512_set1_epi16((int16_t)z[i]);

    uint32_t d = 1;

    for (; d + 31 <= dmax; d += 32) {
      const uint32_t j_low = i - (d + 31);

      __m512i vxj = _mm512_loadu_si512((const void*)(x + j_low));
      __m512i vyj = _mm512_loadu_si512((const void*)(y + j_low));
      __m512i vzj = _mm512_loadu_si512((const void*)(z + j_low));

      __m512i vdx16 = _mm512_sub_epi16(vxi, vxj);
      __m512i vdy16 = _mm512_sub_epi16(vyi, vyj);
      __m512i vdz16 = _mm512_sub_epi16(vzi, vzj);

      // Lower 16 lanes: dt = d+31 .. d+16 (reverse order)
      __m256i dx_lo16 = _mm512_castsi512_si256(vdx16);
      __m256i dy_lo16 = _mm512_castsi512_si256(vdy16);
      __m256i dz_lo16 = _mm512_castsi512_si256(vdz16);

      // Upper 16 lanes: dt = d+15 .. d (reverse order)
      __m256i dx_hi16 = _mm512_extracti64x4_epi64(vdx16, 1);
      __m256i dy_hi16 = _mm512_extracti64x4_epi64(vdy16, 1);
      __m256i dz_hi16 = _mm512_extracti64x4_epi64(vdz16, 1);

      auto widen_store16 = [&](const __m256i& d16, const __m256i& e16, const __m256i& f16, uint32_t* out32) {
        __m128i dlo = _mm256_castsi256_si128(d16);
        __m128i dhi = _mm256_extracti128_si256(d16, 1);
        __m128i elo = _mm256_castsi256_si128(e16);
        __m128i ehi = _mm256_extracti128_si256(e16, 1);
        __m128i flo = _mm256_castsi256_si128(f16);
        __m128i fhi = _mm256_extracti128_si256(f16, 1);

        __m256i d0 = _mm256_cvtepi16_epi32(dlo);
        __m256i d1 = _mm256_cvtepi16_epi32(dhi);
        __m256i e0 = _mm256_cvtepi16_epi32(elo);
        __m256i e1 = _mm256_cvtepi16_epi32(ehi);
        __m256i f0 = _mm256_cvtepi16_epi32(flo);
        __m256i f1 = _mm256_cvtepi16_epi32(fhi);

        __m256i s0 = _mm256_add_epi32(
                       _mm256_add_epi32(_mm256_mullo_epi32(d0,d0), _mm256_mullo_epi32(e0,e0)),
                       _mm256_mullo_epi32(f0,f0));
        __m256i s1 = _mm256_add_epi32(
                       _mm256_add_epi32(_mm256_mullo_epi32(d1,d1), _mm256_mullo_epi32(e1,e1)),
                       _mm256_mullo_epi32(f1,f1));

        _mm256_store_si256((__m256i*)(out32 + 0), s0); // lanes 0..7
        _mm256_store_si256((__m256i*)(out32 + 8), s1); // lanes 8..15
      };

      widen_store16(dx_lo16, dy_lo16, dz_lo16, tmp_sd_lo);
      widen_store16(dx_hi16, dy_hi16, dz_hi16, tmp_sd_hi);

      // Produce dt=d..d+31 in increasing order:
      // dt in [d..d+15]  from tmp_sd_hi reversed (15..0)
      // dt in [d+16..d+31] from tmp_sd_lo reversed (15..0)
      for (uint32_t t = 0; t < 32; ++t) {
        const uint32_t dt = d + t;
        const uint32_t j  = i - dt;

        uint32_t sd = (t < 16) ? tmp_sd_hi[15 - t] : tmp_sd_lo[31 - t];
        const uint16_t code = encode_sd_to_u16(sd);

        out.dist[size_t(i) * M + (dt - 1)] = code;

        if (dt > m_skip) {
          if (out.prev[i] == 0xFFFF ||
              better_candidate(code, dt, j, best_prev_dist, best_prev_d, best_prev_j)) {
            best_prev_dist = code;
            best_prev_d = dt;
            best_prev_j = j;
            out.prev[i] = (uint16_t)j;
          }

          if (out.next[j] == 0xFFFF ||
              better_candidate(code, dt, i, best_next_d[j], best_next_didx[j], out.next[j])) {
            best_next_d[j] = code;
            best_next_didx[j] = dt;
            out.next[j] = (uint16_t)i;
          }
        }
      }
    }

    // Tail scalar
    for (; d <= dmax; ++d) {
      const uint32_t j = i - d;
      const int32_t dx = int32_t(x[i]) - int32_t(x[j]);
      const int32_t dy = int32_t(y[i]) - int32_t(y[j]);
      const int32_t dz = int32_t(z[i]) - int32_t(z[j]);
      const uint32_t sd = uint32_t(dx*dx + dy*dy + dz*dz);
      const uint16_t code = encode_sd_to_u16(sd);

      out.dist[size_t(i) * M + (d - 1)] = code;

      if (d > m_skip) {
        if (out.prev[i] == 0xFFFF ||
            better_candidate(code, d, j, best_prev_dist, best_prev_d, best_prev_j)) {
          best_prev_dist = code;
          best_prev_d = d;
          best_prev_j = j;
          out.prev[i] = (uint16_t)j;
        }

        if (out.next[j] == 0xFFFF ||
            better_candidate(code, d, i, best_next_d[j], best_next_didx[j], out.next[j])) {
          best_next_d[j] = code;
          best_next_didx[j] = d;
          out.next[j] = (uint16_t)i;
        }
      }
    }
  }
}

#endif // AVX512BW

// ============================
// Optional: compile-time "best available"
// (no runtime CPUID; you can add your own dispatch if needed)
// ============================

static inline void kernel_best_available(const uint16_t* x, const uint16_t* y, const uint16_t* z,
                                         uint32_t L, Out out) {
#if defined(__AVX512F__) && defined(__AVX512BW__)
  kernel_avx512bw(x, y, z, L, out);
#elif defined(__AVX2__) || (defined(_MSC_VER) && defined(__AVX2__))
  kernel_avx2(x, y, z, L, out);
#else
  kernel_scalar_ref(x, y, z, L, out);
#endif
}

} // namespace distmx