#include "myutils.h"

static const uint64 _1 = 1;
static const uint64 _2 = 2;
static const uint64 _8 = 8;

/*
 * @brief Converts 2D indices (i, j) for the upper triangular part (with diagonal)
 * of an N x N matrix into a single zero-based index k.
 * The upper triangular part includes all elements where 0 <= i <= j < N.
 * @param i The row index (0-based).
 * @param j The column index (0-based).
 * @param N The dimension of the square matrix (N x N).
 * @return The 1D index k.
 */
uint triangle_ij_to_k(uint i, uint j, uint N)
	{
	asserta(i < N);
	asserta(j < N);
	asserta(i <= j);

    // --- Mathematical Derivation for k ---
    // The number of elements in the first 'i' rows (rows 0 to i-1) is the
    // sum of lengths: N + (N-1) + ... + (N - (i-1)).
    // This is the sum of an arithmetic series: Sum = i * (First + Last) / 2
    // Sum = i * (N + (N - (i-1))) / 2 = i * (2N - i + 1) / 2

    // Use uint64 for intermediate offset calculation to prevent overflow
    // if N is large, though the final k should fit in an uint for typical N.
    uint64 offset = (uint64)i * (_2 * N - i + 1) / _2;

    uint k = (uint)(offset + (j - i));
    return k;
	}

/*
 * @brief Converts a single zero-based index k back into 2D indices (i, j)
 * for the upper triangular part (with diagonal) of an N x N matrix.
 * @param k The 1D index (0-based).
 * @param N The dimension of the square matrix (N x N).
 */
void triangle_k_to_ij(uint k, uint N, uint &i, uint &j)
	{
    uint64 total_elements = (uint64)N * (N + 1) / _2;

	asserta(k < total_elements);
    if (k < 0 || k >= total_elements) {
        throw std::out_of_range("Index k is out of bounds for an N x N upper triangle.");
    }

    // --- Finding Row Index 'i' using Quadratic Formula Approximation ---
    // We are looking for the largest integer 'i' such that the offset S(i) <= k.
    // The formula for the offset S(i) is: S(i) = (i^2 - (2N+1)i + 2k) = 0
    // Solving for 'i' using the quadratic formula: i = (-(2N+1) +/- sqrt((2N+1)^2 - 8k)) / -2
    // Which simplifies to: i = ((2N+1) - sqrt((2N+1)^2 - 8k)) / 2

    // Use uint64 for the intermediate values inside the square root.
    uint64 A = _2 * N + _1;
    uint64 discriminant = A * A - _8*k;

    // Calculate 'i' (the row index)
    // The result of the division is floored by casting to (uint64).
    uint64 i_approx = uint64((A - std::sqrt(discriminant)) / _2);
    
    // We must check if i_approx is the correct floor or if we over-approximated
    // due to floating-point imprecision near the boundary.
    i = (uint)i_approx;
    
    // Calculate the offset for the found row 'i'
    uint64 S_i = (uint64)i * (_2 * N - i + 1) / _2;

    // Refine 'i' if the approximation was slightly off
    if (S_i > k)
		{
        i--;
        S_i = (uint64)i * (_2 * N - i + 1) / _2;
		}
    
    // Calculate the column offset (k - S_i)
    uint col_offset = k - (uint)S_i;

    // The column index j is the row index plus the column offset
    j = i + col_offset;
	}

uint triangle_get_k(uint N)
    {
    return triangle_ij_to_k(N-1, N-1, N);
    }

#if 0
void __cmd_test()
	{
	for (uint N = 100; N < 100000; N += 5000)
		{
		ProgressLog("N=%u\n", N);
		uint PairCount = N + (N*(N-1))/2;
		uint Next1 = 0;
		uint Next2 = 0;
		uint PairIndex = 0;
		for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
			{
			uint idx1, idx2;
			triangle_k_to_ij(PairIndex, N, idx1, idx2);
			uint k = triangle_ij_to_k(idx1, idx2, N);
			asserta(idx1 == Next1);
			asserta(idx2 == Next2);
			asserta(k == PairIndex);
			++Next2;
			if (Next2 == N)
				{
				++Next1;
				Next2 = Next1;
				}
			}
		}
	}
#endif