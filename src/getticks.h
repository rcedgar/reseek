#pragma once

#include <immintrin.h>
#ifndef _MSC_VER
#include <x86intrin.h>
#endif

typedef uint64_t TICKS;

static inline TICKS GetClockTicks() {
    unsigned aux;
    _mm_lfence();
    uint64_t t = __rdtscp(&aux);
    _mm_lfence();
    return t;
}