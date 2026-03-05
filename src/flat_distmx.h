#pragma once

/***
Coordinates are stored as uint16_t "ic".

Distance matrix stores uint32_t squared ic
    Saturates at
    ic_max = 2^16 = 65535
    Angstroms = ic_max/10 = 6553

Squared distance sd(i,j) stored for all 1 < |i-j| <= M
    M pairs (j values) for typical i
    <M pairs close to the ends

Flat matrix layout:
    mx[M*i + j - i - 1] where j = i+1, i+2 ... i+M
    k   = M*(i-1) - 1 + j
***/

static inline uint16_t coord2ic(float x) { return uint16_t((x + 1000)*10 + 0.5); }
static inline float ic2coord(uint16_t ic) { return float(ic/10.0f) - 1000; }

static const uint32_t M = 128;  // band width, M values per i

static inline uint32_t banded_ij_to_k(uint32_t i, uint32_t j)
    {
    if (j < i) std::swap(i, j);
    uint32_t offset = j - i;
    return M*i + offset - 1;
    }

static inline void banded_k_to_ij(uint32_t k, uint32_t& i, uint32_t& j) 
    {
    i = k / M;
    uint32_t offset = (k % M) + 1;
    j = i + offset;
    }