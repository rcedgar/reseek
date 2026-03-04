#include "myutils.h"
#include "dss.h"
#include "pdbchain.h"
#include "flat_chain.h"
#include "distmx_kernel.h"
#include "getticks.h"

const uint M = 100;
const uint skip = 12;

static uint16_t coord2ic(float X) { return uint16_t((X + 1000)*10 + 0.5); }
static int16_t coord2ic_signed(float X) { return int16_t((X + 1000)*10 + 0.5); }
static float ic2coord(uint16_t IC) { return float(IC/10.0f) - 1000; }

static uint16_t GetICDist2(
	uint16_t xic1, uint16_t yic1, uint16_t zic1,
	uint16_t xic2, uint16_t yic2, uint16_t zic2)
	{
	int16_t dx = int(xic1) - int(xic2);
	int16_t dy = int(yic1) - int(yic2);
	int16_t dz = int(zic1) - int(zic2);
	 uint32_t sd = dx*dx + dy*dy + dz*dz;
    return (sd > (uint32_t)DIST_SAT) ? DIST_SAT : (uint16_t)sd;
	}

static TICKS s_TotalTicks;
uint cmp_kernel(vector<uint16_t> &xic, vector<uint16_t> &yic, vector<uint16_t> &zic)
	{
	const bool trace = false;
	uint diffs = 0;
	const uint L = SIZE(xic);
	const uint32_t K = L*M;
	uint16_t* distmx = myalloc(uint16_t, K);
	vector<bool> touched(K);
	TICKS t1 = GetClockTicks();
	dist_fill_avx2(xic.data(), yic.data(), zic.data(), L, M, distmx);
	TICKS t2 = GetClockTicks();
	s_TotalTicks += (t2 - t1);
	for (int i = 0; i < int(L); ++i)
		{
		for (int j = 0; j < i; ++j)
			{
			if (abs(i-j) > M)
				continue;
			uint16_t d2 = GetICDist2(
				xic[i], yic[i], zic[i],
				xic[j], yic[j], zic[j]);
			uint32_t k = dmx_ij_to_k(i, j, M);
			asserta(k < K);
			asserta(!touched[k]);
			touched[k] = true;
			uint16_t d2_kernel = distmx[k];
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

static void GetNNs(const PDBChain &Chain, uint i, uint16_t &men, uint16_t &pen,
	uint M, uint m_skip)
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

static float GetDist(const PDBChain &Chain, uint i, uint j)
	{
	const uint L = Chain.GetSeqLength();
	if (i >= L || j >= L)
		return -1;
	return Chain.GetDist(i, j);
	}

static uint cmp_nn(DSS &D, const uint16_t *next, bool is_next)
	{
	uint diffs = 0;
	const uint L = D.GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		uint nn = (is_next ? D.GetPEN(i) : D.GetMEN(i));
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

void cmd_test_dist_mx()
	{
#if 0
	SimpleTest();
	return;
#else
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	FILE *f = CreateStdioFile(opt(output));
	const uint ChainCount = SIZE(Chains);
	DSS D;
	uint diffs = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Testing %u diffs", diffs);
		const PDBChain &Chain = *Chains[ChainIdx];
		//D.Init(Chain);

		const uint L = Chain.GetSeqLength();

		vector<uint16_t> xic;
		vector<uint16_t> yic;
		vector<uint16_t> zic;
		Chain.GetICsxyz(xic, yic, zic);
		const uint32_t K = L*M;
		uint16_t* distmx = myalloc(uint16_t, K);
		TICKS t1 = GetClockTicks();
		dist_fill_avx2(xic.data(), yic.data(), zic.data(), L, M, distmx);
		TICKS t2 = GetClockTicks();
		s_TotalTicks += (t2 - t1);

		diffs += cmp_kernel(xic, yic, zic);

		//uint16_t* prev = myalloc(uint16_t, L);
		//uint16_t* next = myalloc(uint16_t, L);

		//prev_from_dist(distmx, L, M, skip, prev);
		//next_from_dist(distmx, L, M, skip, next);

		//DSS D;
		//D.Init(Chain);
		//cmp_nn(D, prev, false);
		//cmp_nn(D, next, true);
		}
	ProgressLog("%u chains %.0f ticks %u distance diffs\n",
		ChainCount, double(s_TotalTicks), diffs);
	_chkmem();
#endif
	}