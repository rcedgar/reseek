#pragma once
#include <algorithm>
#include <stdint.h>
#include <string.h>

#ifndef PSORT_OMP_MIN_N
#define PSORT_OMP_MIN_N 65536u
#endif

#ifndef PSORT_BLOCK_SIZE
#define PSORT_BLOCK_SIZE 32768u
#endif

#ifndef PSORT_INSERTION_MAX_BREAKS_DIV
#define PSORT_INSERTION_MAX_BREAKS_DIV 64u
#endif

#ifndef PSORT_INSERTION_MAX_BREAKS_ABS
#define PSORT_INSERTION_MAX_BREAKS_ABS 4096u
#endif

template <class ScoreT, bool Ascending>
struct score_order_less
{
	const ScoreT *scores;

	score_order_less(const ScoreT *s) : scores(s) {}

	inline bool operator()(uint32_t a, uint32_t b) const
	{
		const ScoreT sa = scores[a];
		const ScoreT sb = scores[b];

		if (Ascending)
			{
			if (sa < sb) return true;
			if (sa > sb) return false;
			}
		else
			{
			if (sa > sb) return true;
			if (sa < sb) return false;
			}

		// deterministic tie-break
		return a < b;
	}
};

template <class ScoreT, bool Ascending>
static inline bool score_order_leq(const ScoreT *scores, uint32_t a, uint32_t b)
{
	const ScoreT sa = scores[a];
	const ScoreT sb = scores[b];

	if (Ascending)
		{
		if (sa < sb) return true;
		if (sa > sb) return false;
		}
	else
		{
		if (sa > sb) return true;
		if (sa < sb) return false;
		}

	return a <= b;
}

template <class ScoreT, bool Ascending>
static void insertion_sort_order(const ScoreT *scores, uint32_t n, uint32_t *order)
{
	for (uint32_t i = 1; i < n; ++i)
		{
		const uint32_t x = order[i];
		uint32_t j = i;
		while (j > 0)
			{
			const uint32_t y = order[j - 1];
			if (score_order_leq<ScoreT, Ascending>(scores, y, x))
				break;
			order[j] = y;
			--j;
			}
		order[j] = x;
		}
}

template <class ScoreT, bool Ascending>
static uint32_t count_adjacent_breaks(const ScoreT *scores, uint32_t n, const uint32_t *order)
{
	if (n <= 1)
		return 0;

	uint32_t breaks = 0;

	#if defined(_OPENMP)
	#pragma omp parallel for reduction(+:breaks) schedule(static)
	#endif
	for (int64_t i = 1; i < (int64_t)n; ++i)
		{
		const uint32_t a = order[i - 1];
		const uint32_t b = order[i];
		if (!score_order_leq<ScoreT, Ascending>(scores, a, b))
			++breaks;
		}

	return breaks;
}

template <class ScoreT, bool Ascending>
static void merge_runs(
	const ScoreT *scores,
	const uint32_t *src,
	uint32_t *dst,
	uint32_t left,
	uint32_t mid,
	uint32_t right)
{
	uint32_t i = left;
	uint32_t j = mid;
	uint32_t k = left;

	while (i < mid && j < right)
		{
		const uint32_t a = src[i];
		const uint32_t b = src[j];
		if (score_order_leq<ScoreT, Ascending>(scores, a, b))
			dst[k++] = a, ++i;
		else
			dst[k++] = b, ++j;
		}

	while (i < mid)   dst[k++] = src[i++];
	while (j < right) dst[k++] = src[j++];
}

template <class ScoreT, bool Ascending>
static void parallel_sort_order_impl(const ScoreT *scores, uint32_t n, uint32_t *order)
{
	if (n <= 1)
		return;

	// Count local disorder in the reused permutation.
	// If the old order is still close, insertion sort is often much faster
	// than a full parallel sort because it becomes near-linear.
	const uint32_t breaks = count_adjacent_breaks<ScoreT, Ascending>(scores, n, order);
	if (breaks == 0)
		return;

	const uint32_t break_limit1 = n / PSORT_INSERTION_MAX_BREAKS_DIV;
	const uint32_t break_limit2 = PSORT_INSERTION_MAX_BREAKS_ABS;
	const uint32_t break_limit = break_limit1 < break_limit2 ? break_limit1 : break_limit2;

	if (breaks <= break_limit)
		{
		insertion_sort_order<ScoreT, Ascending>(scores, n, order);
		return;
		}

	// Thread-local scratch reused across calls.
	static thread_local uint32_t *tls_tmp = 0;
	static thread_local uint32_t tls_tmp_cap = 0;

	if (tls_tmp_cap < n)
		{
		if (tls_tmp != 0)
			myfree(tls_tmp);
		tls_tmp = myalloc(uint32_t, n);
		tls_tmp_cap = n;
		}

	uint32_t *tmp = tls_tmp;
	score_order_less<ScoreT, Ascending> less_fn(scores);

	// For smaller n, plain std::sort on the existing order can be faster.
	if (n < PSORT_OMP_MIN_N)
		{
		std::sort(order, order + n, less_fn);
		return;
		}

	// Phase 1: sort blocks independently, in parallel.
	const uint32_t block_size = PSORT_BLOCK_SIZE;
	const uint32_t nblocks = (n + block_size - 1) / block_size;

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (int64_t bi = 0; bi < (int64_t)nblocks; ++bi)
		{
		const uint32_t lo = (uint32_t)bi * block_size;
		uint32_t hi = lo + block_size;
		if (hi > n)
			hi = n;
		std::sort(order + lo, order + hi, less_fn);
		}

	// Phase 2: iterative parallel merges.
	uint32_t *src = order;
	uint32_t *dst = tmp;
	uint32_t width = block_size;

	while (width < n)
		{
		const uint32_t pair_span = width << 1;
		const uint32_t npairs = (n + pair_span - 1) / pair_span;

		#if defined(_OPENMP)
		#pragma omp parallel for schedule(static)
		#endif
		for (int64_t pi = 0; pi < (int64_t)npairs; ++pi)
			{
			const uint32_t left = (uint32_t)pi * pair_span;
			uint32_t mid = left + width;
			uint32_t right = left + pair_span;

			if (mid > n) mid = n;
			if (right > n) right = n;

			if (mid >= right)
				{
				// only one run, copy it
				for (uint32_t k = left; k < right; ++k)
					dst[k] = src[k];
				}
			else
				{
				merge_runs<ScoreT, Ascending>(scores, src, dst, left, mid, right);
				}
			}

		uint32_t *swap = src;
		src = dst;
		dst = swap;
		width <<= 1;
		}

	// If final output ended up in tmp, copy back.
	if (src != order)
		memcpy(order, src, n * sizeof(uint32_t));
}

template <class ScoreT>
static inline void QuickSortOrderT(const ScoreT *scores, uint32_t n, uint32_t *order)
{
	parallel_sort_order_impl<ScoreT, true>(scores, n, order);
}

template <class ScoreT>
static inline void QuickSortOrderDescT(const ScoreT *scores, uint32_t n, uint32_t *order)
{
	parallel_sort_order_impl<ScoreT, false>(scores, n, order);
}

inline void QuickSortOrder_Parallel(const float *scores, uint32_t n, uint32_t *order)
{
	QuickSortOrderT(scores, n, order);
}

inline void QuickSortOrderDesc_Parallel(const float *scores, uint32_t n, uint32_t *order)
{
	QuickSortOrderDescT(scores, n, order);
}

inline void QuickSortOrder_Parallel(const double *scores, uint32_t n, uint32_t *order)
{
	QuickSortOrderT(scores, n, order);
}

inline void QuickSortOrderDesc_Parallel(const double *scores, uint32_t n, uint32_t *order)
{
	QuickSortOrderDescT(scores, n, order);
}