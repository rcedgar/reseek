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

static void compare_kernel(const PDBChain &Chain, const distmx::Out out,
	float t, uint64 &ntest, uint32_t &nfail)
	{
	const uint L = Chain.GetSeqLength();
	for (uint i = 0; i < L; ++i)
		{
		for (uint j = 0; j < i; ++j)
			{
			if (abs(int(i)-int(j)) >= distmx::M)
				continue;
			float d = Chain.GetDist(i, j);
			uint32_t k = distmx::ij_to_k(i, j);
			assert(k < K);
			uint16_t dist_code_kernel = out.dist[k];
			float d_kernel = distmx::decode_code_to_distance(dist_code_kernel);
			++ntest;
			if (fabs(d_kernel - d) > t)
				++nfail;
			}
		}
	}

void cmd_test_dist_mx()
	{
	const bool trace = false;
#if 0
	SimpleTest();
	return;
#endif
#if 1
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	FILE *f = CreateStdioFile(opt(output));
	const uint ChainCount = SIZE(Chains);
	DSS D;
	uint64 ntest_ref = 0;
	uint64 ntest_avx2 = 0;
	uint64 ntest_avx5 = 0;
	uint32 nfail_ref = 0;
	uint32 nfail_avx2 = 0;
	uint32 nfail_avx5 = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Testing nfail_ref=%u", nfail_ref);
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
		compare_kernel(Chain, out, maxdiff, ntest_ref, nfail_ref);

		distmx::kernel_avx2(xic.data(), yic.data(), zic.data(), L, out);
		compare_kernel(Chain, out, maxdiff, ntest_avx2, nfail_avx2);

		distmx::kernel_avx512bw(xic.data(), yic.data(), zic.data(), L, out);
		compare_kernel(Chain, out, maxdiff, ntest_avx5, nfail_avx5);

		free_out(out);
		}
	_chkmem();
	CloseStdioFile(f);
	ProgressLog("%u fail / %s (%.3g%%)\n", 
		nfail_ref, Int64ToStr(ntest_ref), GetPct(nfail_ref, double(ntest_ref)));
	ProgressLog("%u fail / %s (%.3g%%)\n",
		nfail_avx2, Int64ToStr(ntest_avx2), GetPct(nfail_avx2, double(ntest_avx2)));
	ProgressLog("%u fail / %s (%.3g%%)\n",
		nfail_avx5, Int64ToStr(ntest_avx5), GetPct(nfail_avx5, double(ntest_avx5)));
#endif
	}