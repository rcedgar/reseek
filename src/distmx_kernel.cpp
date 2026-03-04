#include "myutils.h"
#include "dss.h"
#include "pdbchain.h"
#include "flat_chain.h"
#include "distmx_kernel.h"
#include "getticks.h"

static uint16_t coord2ic(float X) { return uint16_t((X + 1000)*10 + 0.5); }
static int16_t coord2ic_signed(float X) { return int16_t((X + 1000)*10 + 0.5); }
static float ic2coord(uint16_t IC) { return float(IC/10.0f) - 1000; }

static uint32_t GetICDist2(
	uint16_t xic1, uint16_t yic1, uint16_t zic1,
	uint16_t xic2, uint16_t yic2, uint16_t zic2)
	{
	int16_t dx = int(xic1) - int(xic2);
	int16_t dy = int(yic1) - int(yic2);
	int16_t dz = int(zic1) - int(zic2);
	uint32_t sd = dx*dx + dy*dy + dz*dz;
    return sd;
	}

static TICKS s_DistTicks;
static TICKS s_PENTicks;
static TICKS s_MENTicks;

uint cmp_kernel(vector<uint16_t> &xic, vector<uint16_t> &yic, vector<uint16_t> &zic)
	{
	const bool trace = false;
	uint diffs = 0;
	const uint L = SIZE(xic);
	const uint32_t K = L*M;
	uint32_t* distmx = myalloc(uint32_t, K);
	vector<bool> touched(K);
	TICKS t1 = GetClockTicks();
	dist_fill_avx2(xic.data(), yic.data(), zic.data(), L, distmx);
	TICKS t2 = GetClockTicks();
	s_DistTicks += (t2 - t1);
	for (int i = 0; i < int(L); ++i)
		{
		for (int j = 0; j < i; ++j)
			{
			if (abs(i-j) > M)
				continue;
			uint32_t d2 = GetICDist2(
				xic[i], yic[i], zic[i],
				xic[j], yic[j], zic[j]);
			uint32_t k = dmx_ij_to_k(i, j);
			asserta(k < K);
			asserta(!touched[k]);
			touched[k] = true;
			uint32_t d2_kernel = distmx[k];
			bool is_diff = (d2 != d2_kernel);
			if (is_diff) ++diffs;
			if (trace)
				{
				Log( "%d", i);
				Log( "\t%d", j);
				Log( "\t%d", d2);
				Log( "\t%d", d2_kernel);
				Log( "\n");
				}
			}
		}
	return diffs;
	}

static void SimpleTest()
	{
	const uint L = 400;
	//for (uint i = 0; i < L; ++i)
	//	append(
	//		float(rand()%100),
	//		float(rand()%100),
	//		float(rand()%100));

	//vector<uint16_t> xic;
	//vector<uint16_t> yic;
	//vector<uint16_t> zic;
	//for (uint32_t i = 0; i < L; ++i)
	//	xic.push_back(coord2ic(x[i]));
	//for (uint32_t i = 0; i < L; ++i)
	//	yic.push_back(coord2ic(y[i]));
	//for (uint32_t i = 0; i < L; ++i)
	//	zic.push_back(coord2ic(z[i]));

	//cmp_kernel(xic, yic, zic);
	}

static void GetNNs(const PDBChain &Chain, uint i, 
	uint16_t &men, uint16_t &pen)
	{
	men = UINT16_MAX;
	pen = UINT16_MAX;
	const uint L = Chain.GetSeqLength();
	float mmind = FLT_MAX;
	float pmind = FLT_MAX;
	for (uint j = 0; j < L; ++j)
		{
		if (abs(int(i)-int(j)) > int(M))
			continue;
		if (abs(int(i)-int(j)) <= int(m_skip))
			continue;
		float d = Chain.GetDist(i, j);
		if (j < i && d < mmind)
			{
			men = j;
			mmind = d;
			}
		if (j > i && d < pmind)
			{
			pen = j;
			pmind = d;
			}
		}
	}

static void get_pen(const uint32_t* __restrict distmx, uint L,
	uint16_t* __restrict pen)
	{
	for (uint32_t i = 0; i < L; ++i)
		{
		uint32_t mind2 = UINT32_MAX;
		uint16_t minj = UINT16_MAX;
		for (uint32_t j = i+m_skip; j < min(L, i+M); ++j)
			{
			uint k = j*M + (j - i - 1);
			assert(k == dmx_ij_to_k(i, j));
			if (distmx[k] < mind2)
				{
				mind2 = distmx[k];
				minj = j;
				}
			}
		pen[i] = minj;
		}
	}

static float GetDist(const PDBChain &Chain, uint i, uint j)
	{
	const uint L = Chain.GetSeqLength();
	if (i >= L || j >= L)
		return -1;
	return Chain.GetDist(i, j);
	}

static uint GetNN(const PDBChain &Chain, uint i, bool is_next)
	{
	const uint L = Chain.GetSeqLength();
	uint nn = UINT32_MAX;
	float mindist = FLT_MAX;
	if (is_next)
		{
		for (uint j = i + m_skip; j < L; ++j)
			{
			if (j - i >= M)
				continue;
			float d = Chain.GetDist(i, j);
			if (d < mindist)
				{
				nn = j;
				mindist = d;
				}
			}
		}
	else
		{
		for (int j = 0; j < int(i) - int(m_skip); ++j)
			{
			if (i - j >= M)
				continue;
			float d = Chain.GetDist(i, j);
			if (d < mindist)
				{
				nn = j;
				mindist = d;
				}
			}
		}
	return nn;
	}

static uint cmp_nn(const PDBChain &Chain, const uint16_t *next, bool is_next)
	{
	uint diffs = 0;
	const uint L = Chain.GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		uint nn = GetNN(Chain, i, is_next);
		if (nn == UINT_MAX)
			nn = UINT16_MAX;
		uint nn_kernel = next[i];
		if (nn != nn_kernel)
			{
			Log("%s i=%u  nn=%u  nnk=%u\n", is_next ? "next" : "prev", i, nn, nn_kernel);
			++diffs;
			}
		}
	return diffs;
	}

static void log_indexing(uint L)
	{
	for (uint i = 0; i < L; ++i)
		{
		for (uint j = 0; j < L; ++j)
			{
			if (i == j || abs(int(i) - int(j)) > M)
				continue;
			Log("i=%3u  j=%3u  k=%4u\n", i, j, dmx_ij_to_k(i, j));
			}
		}
	}

static void log_mx(const PDBChain &Chain, const uint32_t *distmx)
	{
	const uint L = Chain.GetSeqLength();
	vector<uint16_t> xic;
	vector<uint16_t> yic;
	vector<uint16_t> zic;
	Chain.GetICsxyz(xic, yic, zic);
	for (uint i = 0; i < L; ++i)
		{
		for (uint j = 0; j < L; ++j)
			{
			if (i == j || abs(int(i) - int(j)) > M)
				continue;
			float d = Chain.GetDist(i, j);

			uint16_t icx_i = xic[i];
			uint16_t icy_i = yic[i];
			uint16_t icz_i = zic[i];

			uint16_t icx_j = xic[j];
			uint16_t icy_j = yic[j];
			uint16_t icz_j = zic[j];

			int dx = int(icx_i) - int(icx_j);
			int dy = int(icy_i) - int(icy_j);
			int dz = int(icz_i) - int(icz_j);

			uint32_t d2_32 = dx*dx + dy*dy + dz*dz;
			uint k = dmx_ij_to_k(i, j);
			Log("i=%3u  j=%3u  k=%4u  d2=%5u mx=%5u  d=%.3g d2=%.3g\n", i, j, k, d2_32, distmx[k], d, d*d);
			}
		}
	}

void cmd_test_dist_mx()
	{
	//log_indexing(150);
	//return;
#if 0
	SimpleTest();
	return;
#else
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	FILE *f = CreateStdioFile(opt(output));
	const uint ChainCount = SIZE(Chains);
	DSS D;
	uint diffs_dist = 0;
	uint diffs_pen = 0;
	uint diffs_men = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Testing %u diffs_pen", diffs_pen);
		const PDBChain &Chain = *Chains[ChainIdx];
		//D.Init(Chain);

		const uint L = Chain.GetSeqLength();

		vector<uint16_t> xic;
		vector<uint16_t> yic;
		vector<uint16_t> zic;
		Chain.GetICsxyz(xic, yic, zic);
		const uint32_t K = L*M;
		uint32_t* distmx = myalloc(uint32_t, K);
		TICKS t1 = GetClockTicks();
		dist_fill_avx2(xic.data(), yic.data(), zic.data(), L, distmx);
		TICKS t2 = GetClockTicks();
		s_DistTicks += (t2 - t1);

		//log_mx(Chain, distmx);

		diffs_dist += cmp_kernel(xic, yic, zic);

		uint16_t* pen = myalloc(uint16_t, L);
		uint16_t* men = myalloc(uint16_t, L);

		TICKS t3 = GetClockTicks();
		get_pen(distmx, L, pen);
		TICKS t4 = GetClockTicks();
		s_PENTicks += (t4 - t3);

		diffs_pen += cmp_nn(Chain, pen, true);

		//next_from_dist(distmx, L, M, skip, next);

		//DSS D;
		//D.Init(Chain);
		//cmp_nn(D, prev, false);
		//cmp_nn(D, next, true);
		}
	ProgressLog("%u chains\n", ChainCount);
	ProgressLog("diffs %u %u %u\n", diffs_dist, diffs_pen, diffs_men);
	ProgressLog("%8.3g dist ticks\n", double(s_DistTicks));
	ProgressLog("%8.3g pen ticks\n", double(s_PENTicks));
	_chkmem();
#endif
	}