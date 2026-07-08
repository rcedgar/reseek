#pragma once

uint triangle_ij_to_k(uint i, uint j, uint N);
void triangle_k_to_ij(uint k, uint N, uint &i, uint &j);
uint triangle_get_K(uint N);

static inline uint triangle_ij_to_k2(uint i, uint j, uint N)
    {
    return triangle_ij_to_k(min(i,j), max(i,j), N);
    }
