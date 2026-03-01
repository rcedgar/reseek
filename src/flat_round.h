#pragma once
#include <cstddef>

// Edit this list to try powers of two, etc.
// Interpreted as *element capacities* (not bytes).
static const std::size_t g_flat_caps[] =
{
    64, 96, 128, 160, 192, 256, 320, 384, 448, 512,
    768, 1024, 1536, 2048, 4096
};
static const std::size_t g_flat_caps_n = sizeof(g_flat_caps) / sizeof(g_flat_caps[0]);

inline std::size_t round_capacity(std::size_t n)
{
    for (std::size_t i = 0; i < g_flat_caps_n; ++i)
        if (n <= g_flat_caps[i])
            return g_flat_caps[i];
    return n; // larger than biggest bucket: no rounding
}