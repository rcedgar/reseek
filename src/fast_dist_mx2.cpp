#include "myutils.h"
#include <immintrin.h>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "fast_dist_mx2.h"
#include "dss.h"
#include "pdbchain.h"

static const uint M = 100;
static const BandIndexLite s_bi(M);

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

static uint16_t CoordToIC(float X) { return uint16_t((X + 1000)*10 + 0.5); }
static float ICToCoord(uint16_t IC) { return float(IC/10.0f) - 1000; }

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

	uint16_t *outIC = myalloc(uint16_t, band_K(L, s_bi));
	banded_distances_avx2_u16_xyz(ICs.data(), L, s_bi, outIC);
	const uint band_size = band_K(L, M);
	vector<bool> touched(band_size);
	uint band_counter = 0;
	for (int i = 0; i < L; ++i)
		{
		for (int j = 0; j < i; ++j)
			{
			if (i-j >= M)
				continue;
			++band_counter;
			float d = GetDist(i, j);
			uint16_t dIC = CoordToIC(d);
			uint32_t k = band_ij_to_k(i, j, L, s_bi);
			assert(k < band_size);
			asserta(!touched[k]);
			touched[k] = true;
			uint16_t dIC2 = outIC[k];
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

void cmd_fast_dist_mx()
	{
#if 1
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
		//uint16_t *outIC = myalloc(uint16_t, L*L);
		uint16_t *outIC = myalloc(uint16_t, band_K(L, s_bi));
		banded_distances_avx2_u16_xyz(ICs.data(), L, s_bi, outIC);
		uint band_counter = 0;
		for (int i = 0; i < Li; ++i)
			{
			for (int j = 1; j < i; ++j)
				{
				if (i-j >= M)
					continue;
				++band_counter;
				float d = Chain.GetDist(uint(i), uint(j));
				uint16_t dIC = Chain.CoordToIC(d);
				uint k = band_ij_to_k(i, j, L, s_bi);
				uint16_t dIC2 = outIC[k];
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
		break;
		}
	_chkmem();
	CloseStdioFile(f);
	}
