#include "myutils.h"
#include "dss.h"
#include "pdbchain.h"
#include "flat_chain.h"
#include "distmx_kernel.h"

static uint16_t coord2ic(float X) { return uint16_t((X + 1000)*10 + 0.5); }
static float ic2coord(uint16_t IC) { return float(IC/10.0f) - 1000; }

static vector<float> x, y, z;
static void append(float X, float Y, float Z)
	{
	x.push_back(X);
	y.push_back(Y);
	z.push_back(Z);
	}

static float GetDist(int i, int j)
	{
	float dx = x[i] - x[j];
	float dy = y[i] - y[j];
	float dz = z[i] - z[j];
	return sqrt(dx*dx + dy*dy + dz*dz);
	}

static void alloc_out(distmx::Out &out, uint32_t L)
	{
	out.dist = myalloc(uint16_t, distmx::dist_buffer_length(L));
	out.prev = myalloc(uint16_t, L);
	out.next = myalloc(uint16_t, L);
	}

static void free_out(distmx::Out &out)
	{
	myfree(out.dist);
	myfree(out.prev);
	myfree(out.next);
	}

static void SimpleTest()
	{
	const uint L = 400;
	for (uint i = 0; i < L; ++i)
		append(
			float(rand()%100),
			float(rand()%100),
			float(rand()%100));

	vector<uint16_t> xic;
	vector<uint16_t> yic;
	vector<uint16_t> zic;
	for (uint32_t i = 0; i < L; ++i)
		xic.push_back(coord2ic(x[i]));
	for (uint32_t i = 0; i < L; ++i)
		yic.push_back(coord2ic(y[i]));
	for (uint32_t i = 0; i < L; ++i)
		zic.push_back(coord2ic(z[i]));

	const uint32_t K = distmx::dist_buffer_length(L);
	vector<bool> touched(K);
	distmx::Out out;
	alloc_out(out, L);
	distmx::kernel_scalar_ref(xic.data(), yic.data(), zic.data(), L, out);
	for (int i = 0; i < L; ++i)
		{
		for (int j = 0; j < i; ++j)
			{
			if (abs(i-j) >= distmx::M)
				continue;
			float d = GetDist(i, j);
			uint32_t k = distmx::ij_to_k(i, j);
			assert(k < K);
			asserta(!touched[k]);
			touched[k] = true;
			uint16_t dist_code_kernel = out.dist[k];
			float d_kernel = distmx::decode_code_to_distance(dist_code_kernel);

			Log( "%d", i);
			Log( "\t%d", j);
			Log( "\t%.3g", d);
			Log( "\t%.3g", d_kernel);
			Log( "\t%.3g", d_kernel-d);
			Log( "\n");
			}
		}
	}

static void GetNNs(const PDBChain &Chain, uint i, uint16_t &men, uint16_t &pen)
	{
	men = UINT16_MAX;
	pen = UINT16_MAX;
	const uint L = Chain.GetSeqLength();
	float mmind = FLT_MAX;
	float pmind = FLT_MAX;
	for (uint j = 0; j < L; ++j)
		{
		if (abs(int(i)-int(j)) > distmx::M)
			continue;
		if (abs(int(i)-int(j)) <= distmx::m_skip)
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

static void compare_kernel(const PDBChain &Chain, const distmx::Out out,
	float t, uint64 &ndtest, uint32_t &ndfail, uint32_t &nntest, uint32_t &nnfail)
	{
	const uint L = Chain.GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		uint16_t men, pen;
		GetNNs(Chain, i, men, pen);
		uint16_t men_kernel = out.prev[i];
		uint16_t pen_kernel = out.next[i];
		++nntest;
		if (men_kernel != men) 
			{
			++nnfail;
			Log("[%4u]  %4u(%7.3g) %4u(%7.3g) MEN diff\n",
				i,
				men, GetDist(Chain, i, men),
				men_kernel, GetDist(Chain, i, men_kernel));
			}
		if (pen_kernel != pen)
			{
			++nnfail;
			Log("[%4u]  %4u(%7.3g) %4u(%7.3g) PEN diff\n",
				i,
				pen, GetDist(Chain, i, pen),
				pen_kernel, GetDist(Chain, i, pen_kernel));
			}

		for (uint j = 0; j < i; ++j)
			{
			if (abs(int(i)-int(j)) >= distmx::M)
				continue;
			float d = Chain.GetDist(i, j);
			uint32_t k = distmx::ij_to_k(i, j);
			assert(k < K);
			uint16_t dist_code_kernel = out.dist[k];
			float d_kernel = distmx::decode_code_to_distance(dist_code_kernel);
			++ndtest;
			if (fabs(d_kernel - d) > t)
				++ndfail;
			}
		}
	}

void cmd_test_dist_mx()
	{
	const bool trace = false;
#if 0
	SimpleTest();
	return;
#else
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	FILE *f = CreateStdioFile(opt(output));
	const uint ChainCount = SIZE(Chains);
	DSS D;
	uint64 ndtest_ref = 0;
	uint64 ndtest_avx2 = 0;
	uint64 ndtest_avx5 = 0;
	uint32 ndfail_ref = 0;
	uint32 ndfail_avx2 = 0;
	uint32 ndfail_avx5 = 0;
	uint32 nntest_ref = 0;
	uint32 nnfail_ref = 0;
	uint32 nntest_avx2 = 0;
	uint32 nnfail_avx2 = 0;
	uint32 nntest_avx5 = 0;
	uint32 nnfail_avx5 = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Testing nnfail_ref=%u", nnfail_ref);
		const PDBChain &Chain = *Chains[ChainIdx];
		D.Init(Chain);

		const uint L = Chain.GetSeqLength();
		const uint32_t K = distmx::dist_buffer_length(L);
		const int iL = L;
		distmx::Out out;
		alloc_out(out, L);

		vector<uint16_t> xic;
		vector<uint16_t> yic;
		vector<uint16_t> zic;
		Chain.GetICsxyz(xic, yic, zic);

		const float maxdiff = 0.01;

		distmx::kernel_scalar_ref(xic.data(), yic.data(), zic.data(), L, out);
		compare_kernel(Chain, out, maxdiff, ndtest_ref, ndfail_ref, nntest_ref, nnfail_ref);

		distmx::kernel_avx2(xic.data(), yic.data(), zic.data(), L, out);
		compare_kernel(Chain, out, maxdiff, ndtest_avx2, ndfail_avx2, nntest_avx2, nnfail_avx2);

		distmx::kernel_avx512bw(xic.data(), yic.data(), zic.data(), L, out);
		compare_kernel(Chain, out, maxdiff, ndtest_avx5, ndfail_avx5, nntest_avx5, nnfail_avx5);

		free_out(out);
		}
	_chkmem();
	CloseStdioFile(f);
	ProgressLog("%u dist fail / %s (%.3g%%)\n", 
		ndfail_ref, Int64ToStr(ndtest_ref), GetPct(ndfail_ref, double(ndtest_ref)));
	ProgressLog("%u dist fail / %s (%.3g%%)\n",
		ndfail_avx2, Int64ToStr(ndtest_avx2), GetPct(ndfail_avx2, double(ndtest_avx2)));
	ProgressLog("%u dist fail / %s (%.3g%%)\n",
		ndfail_avx5, Int64ToStr(ndtest_avx5), GetPct(ndfail_avx5, double(ndtest_avx5)));

	ProgressLog("%u nn fail / %s (%.3g%%)\n", 
		nnfail_ref, Int64ToStr(nntest_ref), GetPct(nnfail_ref, double(nntest_ref)));
	ProgressLog("%u nn fail / %s (%.3g%%)\n",
		nnfail_avx2, Int64ToStr(nntest_avx2), GetPct(nnfail_avx2, double(nntest_avx2)));
	ProgressLog("%u nn fail / %s (%.3g%%)\n",
		nnfail_avx5, Int64ToStr(nntest_avx5), GetPct(nnfail_avx5, double(nntest_avx5)));
#endif
	}