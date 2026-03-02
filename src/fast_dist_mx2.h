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
