#include "myutils.h"
#include "fast_dist_mx.h"
#include "dss.h"

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

static void SimpleTest()
	{
	const uint L = 400;
	for (uint i = 0; i < L; ++i)
		append(
			float(rand()%100),
			float(rand()%100),
			float(rand()%100));

	vector<uint16_t> ICs;
	ICs.reserve(3*L);
	for (uint32_t i = 0; i < L; ++i)
		ICs.push_back(CoordToIC(x[i]));
	for (uint32_t i = 0; i < L; ++i)
		ICs.push_back(CoordToIC(y[i]));
	for (uint32_t i = 0; i < L; ++i)
		ICs.push_back(CoordToIC(z[i]));

	uint16_t *outIC = myalloc(uint16_t, band_K(L));
	uint16_t *men = myalloc(uint16_t, L);
	uint16_t *pen = myalloc(uint16_t, L);
	band2_distances_avx2_u16_xyz_v5(ICs.data(), L, outIC, men, pen);
	const uint band_size = band_K(L);
	vector<bool> touched(band_size);
	uint band_counter = 0;
	for (int i = 0; i < L; ++i)
		{
		for (int j = 0; j < i; ++j)
			{
			if (i-j >= dist_mx_band_width)
				continue;
			++band_counter;
			float d = GetDist(i, j);
			uint16_t dIC = CoordToIC(d);
			uint32_t k = band_ij_to_k(i, j);
			assert(k < band_size);
			asserta(!touched[k]);
			touched[k] = true;
			uint16_t dIC2 = outIC[k];
			uint16_t dIC3 = band2_get_checked(outIC, i, j);
			asserta(dIC3 == dIC2);
			int diff = int(dIC) - int(dIC2);
			Log( "%d", i);
			Log( "\t%d", j);
			Log( "\t%u", dIC);
			Log( "\t%u", dIC2);
			if (diff != 0) Log( "\t%d", diff);
			Log( "\n");
			}
		}
	uint not_touched = 0;
	for (uint i = 0; i < band_size; ++i)
		if (!touched[i])
			++not_touched;
	Log("K=%u, count=%u, not=%u\n", band_size, band_counter, not_touched);
	}

#if 0
void cmd_fast_dist_mx()
	{
#if 0
	SimpleTest();
	return;
#endif
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	FILE *f = CreateStdioFile(opt(output));
	const uint ChainCount = SIZE(Chains);
	DSS D;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *Chains[ChainIdx];
		D.Init(Chain);
		vector<uint16_t> ICs;
		Chain.GetICsxyz(ICs);

		const uint L = Chain.GetSeqLength();
		const int Li = L;
		uint32_t K = band_K(L);
		uint16_t *outIC = myalloc(uint16_t, K);
		uint16_t *men = myalloc(uint16_t, L);
		uint16_t *pen = myalloc(uint16_t, L);
		band2_distances_avx2_u16_xyz_v5(ICs.data(), L, outIC, men, pen);
		uint band_counter = 0;
		for (int i = 0; i < Li; ++i)
			{
			for (int j = 1; j < i; ++j)
				{
				if (i-j >= dist_mx_band_width)
					continue;
				++band_counter;
				float d = Chain.GetDist(uint(i), uint(j));
				uint16_t dIC = Chain.CoordToIC(d);
				uint k = band_ij_to_k(i, j);
				uint16_t dIC2 = outIC[k];
				uint16_t dIC3 = band2_get_checked(outIC, i, j);
				asserta(dIC3 == dIC2);
				int diff = int(dIC2) - int(dIC);
				fprintf(f, "%d", i);
				fprintf(f, "\t%d", j);
				fprintf(f, "\t%u", dIC);
				fprintf(f, "\t%u", dIC2);
				if (diff != 0)
					fprintf(f, "\t%d", diff);
				fprintf(f, "\n");
				}
			}
		_chkmem();
		}
	CloseStdioFile(f);
	}
#endif // 0

void cmd_fast_dist_mx()
	{
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	const uint ChainCount = SIZE(Chains);
	const uint M = 100;
	const uint m = 12;
	uint64_t N = 0;
	uint64_t n = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "n=%s", Int64ToStr(n));
		const PDBChain &Chain = *Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		const int iL = int(L);
		for (int i = 0; i < iL ; ++i)
			{
			for (int j = 1; j < i; ++j)
				{
				int dij = abs(i - j);
				if (dij < m || dij > M)
					continue;
				float d2 = Chain.GetDist2(i, j);
				uint32_t ic = uint32_t(d2 + 0.5);
				++N;
				if (ic >= UINT16_MAX)
					++n;
				}
			}
		}
	ProgressLog("N=%s", Int64ToStr(N));
	ProgressLog(" n=%s", Int64ToStr(n));
	ProgressLog("\n");
	}
