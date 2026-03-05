#include "myutils.h"
#include "pdbchain.h"
#include "flat_distmx.h"
#include "getticks.h"

void test_indexing()
	{
	const int L = 250;
	const uint K = 250*M;
	vector<bool> touched_plus(K);
	for (int i = 0; i < L; ++i)
		{
		for (int j = i+1; j < L; ++j)
			{
			if (abs(i-j) > M)
				continue;
			uint32_t k = banded_ij_to_k(i, j);
			asserta(k < K);

			asserta(!touched_plus[k]);
			touched_plus[k] = true;

			uint32_t i2, j2;
			banded_k_to_ij(k, i2, j2);
			if (i2 != i)
				Die("i2=%u i=%u j=%u", i2, i, j);
			if (j2 != j)
				Die("j2=%u i=%u j=%u", j2, i, j);
			}
		}

	vector<bool> touched_minus(K);
	for (int i = 0; i < L; ++i)
		{
		for (int j = i+1; j < L; ++j)
			{
			if (abs(i-j) > M)
				continue;
			uint32_t k = banded_ij_to_k(i, j);
			asserta(k < K);

			asserta(!touched_minus[k]);
			touched_minus[k] = true;
			uint32_t i2, j2;
			banded_k_to_ij(k, i2, j2);
			if (i2 != i)
				Die("i2=%u i=%u j=%u", i2, i, j);
			if (j2 != j)
				Die("j2=%u i=%u j=%u", j2, i, j);
			}
		}
	ProgressLog("test_indexing OK\n");
	}

#if TRACE
static const uint trace_i = 0;
static const uint trace_j = 1;
#endif

static inline void fill_flat_distmx(
	const uint16_t *__restrict xyz,
	uint32_t L,
	uint32_t *__restrict sdmx)
	{
	uint i3 = 0;
	for (uint32_t i = 0; i < L; ++i)
		{
		int32_t icx_i = xyz[i3++];
		int32_t icy_i = xyz[i3++];
		int32_t icz_i = xyz[i3++];
		uint32_t k = i*M;
		const uint32_t jend = min(i+M, L-1);
		for (uint32_t j = i + 1; j <= jend; ++j)
			{
			int32_t dicx = icx_i - xyz[3*j];
			int32_t dicy = icy_i - xyz[3*j+1];
			int32_t dicz = icz_i - xyz[3*j+2];
			assert(k == banded_ij_to_k(i, j));
			uint32_t sd = dicx*dicx + dicy*dicy + dicz*dicz;
			sdmx[k++] = sd;
			}
		}
	}

static inline void fill_pen(uint32_t *__restrict sdmx,
	uint32_t L, uint32_t m, uint16_t *pen)
	{
	for (uint32_t i = 0; i < L; ++i)
		{
		const uint32_t jend = min(i+M, L-1);
		uint16_t pen_i = UINT16_MAX;
		uint32_t min_sd = UINT32_MAX;
		uint32_t k = i*M + m - 1;
		for (uint32_t j = i + m; j <= jend; ++j)
			{
			assert(k == banded_ij_to_k(i, j));
			uint32_t sd = sdmx[k++];
			if (sd < min_sd)
				{
				min_sd = sd;
				pen_i = j;
				}
			}
		pen[i] = pen_i;
		}
	}

static inline void fill_men(uint32_t *__restrict sdmx,
	uint32_t L, uint32_t m, uint16_t *men)
	{
	for (int i = 0; i < int(L); ++i)
		{
		const uint32_t jstart = min(i+M, L-1);
		uint16_t men_i = UINT16_MAX;
		uint32_t min_sd = UINT32_MAX;
		for (int j = max(0,i-int(M)); j <= i-int(m); ++j)
			{
			uint32_t k = banded_ij_to_k(i, j);
			uint32_t sd = sdmx[k];
			if (sd < min_sd)
				{
				min_sd = sd;
				men_i = j;
				}
			}
		men[i] = men_i;
		}
	}

static uint compare_fill(const PDBChain &Chain, const uint32_t *sdmx)
	{
	const uint L = Chain.GetSeqLength();
	vector<uint16_t> x;
	vector<uint16_t> y;
	vector<uint16_t> z;
	Chain.GetICsxyz(x, y, z);

	vector<uint16_t> ICs;
	Chain.GetICs(ICs);

	for (uint i = 0; i < L; ++i)
		{
		asserta(x[i] == ICs[3*i]);
		asserta(y[i] == ICs[3*i+1]);
		asserta(z[i] == ICs[3*i+2]);
		}

	uint diffs = 0;
	for (int i = 0; i < int(L); ++i)
		{
		for (int j = i+1; j < int(L); ++j)
			{
			if (i==j || abs(i-j) > M)
				continue;
			int dx = int(x[i]) - int(x[j]);
			int dy = int(y[i]) - int(y[j]);
			int dz = int(z[i]) - int(z[j]);
			uint32_t sd = dx*dx + dy*dy + dz*dz;
			uint k = banded_ij_to_k(i, j);
			uint32_t sd2 = sdmx[k];
			if (sd2 != sd)
				++diffs;
			}
		}
	return diffs;
	}

static uint compare_pen(const PDBChain &Chain, uint m, const uint16_t *pen)
	{
	const uint L = Chain.GetSeqLength();
	vector<uint16_t> x;
	vector<uint16_t> y;
	vector<uint16_t> z;
	Chain.GetICsxyz(x, y, z);

	uint diffs = 0;
	for (int i = 0; i < int(L); ++i)
		{
		const int jend = min(i+int(M), int(L)-1);
		float MinDist = FLT_MAX;
		uint16_t pen_i = UINT16_MAX;
		for (int j = i + int(m); j <= jend; ++j)
			{
			if (i==j || abs(i-j) < int(m) || abs(i-j) > M)
				continue;
			float d = Chain.GetDist(i, j);
			if (d < MinDist)
				{
				MinDist = d;
				pen_i = j;
				}
			}
		if (pen_i != pen[i])
			++diffs;
		}
	return diffs;
	}

static uint compare_men(const PDBChain &Chain, uint m, const uint16_t *men)
	{
	const uint L = Chain.GetSeqLength();
	vector<uint16_t> x;
	vector<uint16_t> y;
	vector<uint16_t> z;
	Chain.GetICsxyz(x, y, z);

	uint diffs = 0;
	for (int i = 0; i < int(L); ++i)
		{
		const int jend = min(i+int(M), int(L)-1);
		float MinDist = FLT_MAX;
		uint16_t men_i = UINT16_MAX;
		for (int j = 0; j < i; ++j)
			{
			if (i==j || abs(i-j) < int(m) || abs(i-j) > M)
				continue;
			float d = Chain.GetDist(i, j);
			if (d < MinDist)
				{
				MinDist = d;
				men_i = j;
				}
			}
		if (men_i != men[i])
			++diffs;
		}
	return diffs;
	}

static void test_sd(const vector<PDBChain *> &Chains)
	{
	const uint ChainCount = SIZE(Chains);
	vector<uint16_t> ICs;
	uint total_diffs = 0;
	TICKS total_ticks = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "working diffs %u", total_diffs);
		const PDBChain &Chain = *Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		const uint K = L*M;
		Chain.GetICs(ICs);
		const uint16_t *xyz = ICs.data();
		uint32_t *sdmx = myalloc(uint32_t, K);
		TICKS t1 = GetClockTicks();
		fill_flat_distmx(xyz, L, sdmx);
		TICKS t2 = GetClockTicks();
		total_ticks += t2 - t1;
		uint diffs = compare_fill(Chain, sdmx);
		total_diffs += diffs;
		}
	ProgressLog("%.3g ticks, %u diffs sd\n", double(total_ticks), total_diffs);
	}

static void test_pen(const vector<PDBChain *> &Chains, uint m)
	{
	const uint ChainCount = SIZE(Chains);
	vector<uint16_t> ICs;
	uint total_diffs = 0;
	TICKS total_ticks = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "working diffs %u", total_diffs);
		const PDBChain &Chain = *Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		const uint K = L*M;
		Chain.GetICs(ICs);
		const uint16_t *xyz = ICs.data();
		uint32_t *sdmx = myalloc(uint32_t, K);
		fill_flat_distmx(xyz, L, sdmx);

		uint16_t *pen = myalloc(uint16_t, L);

		TICKS t1 = GetClockTicks();
		fill_pen(sdmx, L, m, pen);
		TICKS t2 = GetClockTicks();

		total_ticks += t2 - t1;
		uint diffs = compare_pen(Chain, m, pen);
		total_diffs += diffs;
		}
	ProgressLog("%.3g ticks, %u diffs pen\n", double(total_ticks), total_diffs);
	}

static void test_men(const vector<PDBChain *> &Chains, uint m)
	{
	const uint ChainCount = SIZE(Chains);
	vector<uint16_t> ICs;
	uint total_diffs = 0;
	TICKS total_ticks = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "working diffs %u", total_diffs);
		const PDBChain &Chain = *Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		const uint K = L*M;
		Chain.GetICs(ICs);
		const uint16_t *xyz = ICs.data();
		uint32_t *sdmx = myalloc(uint32_t, K);
		fill_flat_distmx(xyz, L, sdmx);

		uint16_t *men = myalloc(uint16_t, L);

		TICKS t1 = GetClockTicks();
		fill_men(sdmx, L, m, men);
		TICKS t2 = GetClockTicks();

		total_ticks += t2 - t1;
		uint diffs = compare_men(Chain, m, men);
		total_diffs += diffs;
		}
	ProgressLog("%.3g ticks, %u diffs men\n", double(total_ticks), total_diffs);
	}

void cmd_test_flat_distmx()
	{
	test_indexing();
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	test_sd(Chains);
	test_pen(Chains, 16);
	test_men(Chains, 16);
	}