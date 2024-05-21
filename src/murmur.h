#pragma once

static inline uint32_t murmur_32_scramble(uint32_t k)
    {
    k *= 0xcc9e2d51;
    k = (k << 15) | (k >> 17);
    k *= 0x1b873593;
    return k;
    }

uint32_t murmur3_32(const byte* key, size_t len,
  uint32_t seed = 0x9747b28c);
