#include "myutils.h"
#include <immintrin.h>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "fast_dist_mx2.h"
#include "dss.h"
#include "pdbchain.h"

static inline void banded_distances_avx2_u16_xyz(
	const uint16_t* __restrict xyz,   // length 3*N: [x...][y...][z...], each uint16 IC
	uint32_t N,
	const BandIndexLite& bi,          // constructed once for your fixed M (e.g. M=100)
	uint16_t* __restrict outIC        // length band_K(N, bi)
){
	assert(xyz && outIC);
	if (N < 2 || bi.w == 0) return;

	const uint16_t* __restrict x = xyz;
	const uint16_t* __restrict y = xyz + N;
	const uint16_t* __restrict z = xyz + 2u*N;

	const __m256i i10000 = _mm256_set1_epi32(10000);
	const __m256i izero  = _mm256_setzero_si256();
	const __m256i i65535 = _mm256_set1_epi32(65535);

	for (uint32_t i = 1; i < N; ++i) {
		uint32_t len = std::min(i, bi.w);
		if (len == 0) continue;

		uint32_t start = i - len;                 // max(0, i-w)
		uint32_t base  = band_base(i, bi);        // row base in flattened output
		uint32_t t = 0;

		const __m256i xi = _mm256_set1_epi32((int32_t)x[i]);
		const __m256i yi = _mm256_set1_epi32((int32_t)y[i]);
		const __m256i zi = _mm256_set1_epi32((int32_t)z[i]);

		for (; t + 7 < len; t += 8) {
			uint32_t j = start + t;

			__m128i x8u16 = _mm_loadu_si128((const __m128i*)(x + j));
			__m128i y8u16 = _mm_loadu_si128((const __m128i*)(y + j));
			__m128i z8u16 = _mm_loadu_si128((const __m128i*)(z + j));

			__m256i xj = _mm256_cvtepu16_epi32(x8u16);
			__m256i yj = _mm256_cvtepu16_epi32(y8u16);
			__m256i zj = _mm256_cvtepu16_epi32(z8u16);

			__m256i dx = _mm256_sub_epi32(xj, xi);
			__m256i dy = _mm256_sub_epi32(yj, yi);
			__m256i dz = _mm256_sub_epi32(zj, zi);

			__m256i d2 = _mm256_mullo_epi32(dx, dx);
			d2 = _mm256_add_epi32(d2, _mm256_mullo_epi32(dy, dy));
			d2 = _mm256_add_epi32(d2, _mm256_mullo_epi32(dz, dz));

			__m256 fs = _mm256_sqrt_ps(_mm256_cvtepi32_ps(d2));

			// deterministic round-to-nearest: round then truncate
			fs = _mm256_round_ps(fs, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
			__m256i s_int = _mm256_cvttps_epi32(fs);

			__m256i dic = _mm256_add_epi32(s_int, i10000);
			dic = _mm256_max_epi32(dic, izero);
			dic = _mm256_min_epi32(dic, i65535);

			// pack 8x i32 -> 8x u16 contiguously (correct)
			__m128i dic_lo = _mm256_castsi256_si128(dic);
			__m128i dic_hi = _mm256_extracti128_si256(dic, 1);
			__m128i out8u16 = _mm_packus_epi32(dic_lo, dic_hi);

			_mm_storeu_si128((__m128i*)(outIC + base + t), out8u16);
		}

		for (; t < len; ++t) {
			uint32_t j = start + t;

			int32_t dx = (int32_t)x[j] - (int32_t)x[i];
			int32_t dy = (int32_t)y[j] - (int32_t)y[i];
			int32_t dz = (int32_t)z[j] - (int32_t)z[i];
			int32_t d2 = dx*dx + dy*dy + dz*dz;

			float s = std::sqrt((float)d2);                 // deci-units
			uint32_t si = (uint32_t)std::floor(s + 0.5f);   // round-to-nearest (s>=0)
			uint32_t dic = si + 10000u;
			if (dic > 65535u) dic = 65535u;
			outIC[base + t] = (uint16_t)dic;
		}
	}
}

static BandIndexLite s_bi(1000);

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
static const uint M = 100;

static void SimpleTest()
	{
	const uint L = 50;
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
	for (int i = 0; i < L; ++i)
		{
		for (int j = 1; j < i; ++j)
			{
			if (i-j >= M)
				break;
			float d = GetDist(i, j);
			uint16_t dIC = CoordToIC(d);
			uint32_t k = band_ij_to_k(i, j, L, s_bi);
			uint16_t dIC2 = outIC[k];
			printf( "%d", i);
			printf( "\t%d", j);
			printf( "\t%u", dIC);
			printf( "\t%u", dIC2);
			printf( "\n");
			}
		}
	}

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
		uint16_t *outIC = myalloc(uint16_t, L*L);
		banded_distances_avx2_u16_xyz(ICs.data(), L, s_bi, outIC);

		for (int i = 0; i < Li; ++i)
			{
			for (int j = 1; j < i; ++j)
				{
				if (i-j >= M)
					break;
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
	CloseStdioFile(f);
	}
