#pragma once

#include <vector>
#include <cstdint>
#include <cassert>

template <typename T, bool GetMin>
std::vector<uint32_t> get_distinct_window_extrema(const T* v, uint32_t n, uint32_t w)
	{
	std::vector<uint32_t> idxs;

	if (v == nullptr || n == 0 || w == 0)
		return idxs;

	if (w > n)
		w = n;

	// Mark positions that are selected as the extrema of at least one window.
	std::vector<uint8_t> is_extreme(n, 0);

	const uint32_t nwin = n - w + 1;
	for (uint32_t winlo = 0; winlo < nwin; ++winlo)
		{
		uint32_t besti = winlo;
		T best = v[winlo];

		const uint32_t winhi = winlo + w;
		for (uint32_t i = winlo + 1; i < winhi; ++i)
			{
			const T x = v[i];
			if constexpr (GetMin)
				{
				if (x < best)
					{
					best = x;
					besti = i;
					}
				}
			else
				{
				if (x > best)
					{
					best = x;
					besti = i;
					}
				}
			}

		is_extreme[besti] = 1;
		}

	for (uint32_t i = 0; i < n; ++i)
		{
		if (is_extreme[i])
			idxs.push_back(i);
		}

	return idxs;
	}
