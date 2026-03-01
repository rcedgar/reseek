#include <cstdint>
#include <vector>
#include <algorithm>
#include <immintrin.h>
#include <cmath>
#include <cassert>
#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"

// https://chatgpt.com/c/69a323b2-3a08-832c-8180-ee5e11e53c61

static const int M = 100;

static uint16_t CoordToIC(float X) { return uint16_t((X + 1000)*10 + 0.5); }
static float ICToCoord(uint16_t IC) { return float(IC/10.0f) - 1000; }

struct BandIndex {
	uint32_t N = 0;
	uint32_t M = 0;                 // abs(i-j) < M, so lower band width is (M-1)
	std::vector<uint32_t> base;     // size N+1, base[0]=0, base[N]=K

	uint32_t K() const { return base.empty() ? 0u : base[N]; }
};

// base[i] = total entries in rows [0..i-1]
static inline BandIndex band_build_index(uint32_t N, uint32_t M) {
	BandIndex idx;
	idx.N = N;
	idx.M = M;
	idx.base.resize(N + 1);
	idx.base[0] = 0;
	uint32_t w = (M > 0 ? M - 1 : 0);          // lower-band max distance in index units
	for (uint32_t i = 0; i < N; ++i) {
		uint32_t len = std::min(i, w);
		idx.base[i + 1] = idx.base[i] + len;
	}
	return idx;
}

static inline uint32_t band_ij_to_k(uint32_t i, uint32_t j, const BandIndex& idx) {
    assert(i < idx.N && j < idx.N);
    assert(j < i);                              // lower only
    uint32_t w = (idx.M > 0 ? idx.M - 1 : 0);
    assert(i - j <= w);                         // band condition (i-j < M)

    uint32_t len   = std::min(i, w);
    uint32_t start = i - len;                   // max(0, i-w)
    return idx.base[i] + (j - start);
}

// This is not in the hot path; do a binary search in base[] (N <= 500, trivial).
static inline void band_k_to_ij(uint32_t k, const BandIndex& idx, uint32_t &i, uint32_t &j) {
    assert(k < idx.K());
    // Find i such that base[i] <= k < base[i+1]
    auto it = std::upper_bound(idx.base.begin(), idx.base.end(), k);
    i = uint32_t((it - idx.base.begin()) - 1);

    uint32_t w   = (idx.M > 0 ? idx.M - 1 : 0);
    uint32_t len = std::min(i, w);
    uint32_t start = i - len;
    uint32_t t = k - idx.base[i];
    j = start + t;

    assert(j < i);
    assert(i - j <= w);
}

// Compute banded lower-triangle distances:
// store pairs (i,j) with j<i and (i-j)<M
// Output distances are IC format: round(dist*10) + 10000 (since dist>=0).
static inline void banded_distances_avx2_u16(
    const uint16_t* __restrict coordsIC,   // length 3*N: x0,y0,z0,x1,y1,z1,...
    uint32_t N,
    uint32_t M,
    uint16_t* __restrict outIC,            // length idx.K() (you may allocate bigger)
    const BandIndex& idx
) {
    assert(idx.N == N && idx.M == M);
    assert(outIC && coordsIC);

    // AoS -> SoA
    std::vector<uint16_t> xs(N), ys(N), zs(N);
    for (uint32_t i = 0; i < N; ++i) {
        xs[i] = coordsIC[3*i + 0];
        ys[i] = coordsIC[3*i + 1];
        zs[i] = coordsIC[3*i + 2];
    }

    const uint32_t w = (M > 0 ? M - 1 : 0);
    const __m256i i10000 = _mm256_set1_epi32(10000);
    const __m256i izero  = _mm256_setzero_si256();
    const __m256i i65535 = _mm256_set1_epi32(65535);

    for (uint32_t i = 1; i < N; ++i) {
        uint32_t len = std::min(i, w);
        if (len == 0) continue;

        uint32_t start = i - len;     // max(0, i-w)
        uint32_t base  = idx.base[i]; // output row base
        uint32_t t = 0;

        const __m256i xi = _mm256_set1_epi32((int32_t)xs[i]);
        const __m256i yi = _mm256_set1_epi32((int32_t)ys[i]);
        const __m256i zi = _mm256_set1_epi32((int32_t)zs[i]);

        for (; t + 7 < len; t += 8) {
            uint32_t j = start + t;

            __m128i x8u16 = _mm_loadu_si128((const __m128i*)(xs.data() + j));
            __m128i y8u16 = _mm_loadu_si128((const __m128i*)(ys.data() + j));
            __m128i z8u16 = _mm_loadu_si128((const __m128i*)(zs.data() + j));

            __m256i xj = _mm256_cvtepu16_epi32(x8u16);
            __m256i yj = _mm256_cvtepu16_epi32(y8u16);
            __m256i zj = _mm256_cvtepu16_epi32(z8u16);

            __m256i dx = _mm256_sub_epi32(xj, xi);
            __m256i dy = _mm256_sub_epi32(yj, yi);
            __m256i dz = _mm256_sub_epi32(zj, zi);

            __m256i d2 = _mm256_mullo_epi32(dx, dx);
            d2 = _mm256_add_epi32(d2, _mm256_mullo_epi32(dy, dy));
            d2 = _mm256_add_epi32(d2, _mm256_mullo_epi32(dz, dz));

            // sqrt(d2) in float (d2 is in deci-units^2)
            __m256 fs = _mm256_sqrt_ps(_mm256_cvtepi32_ps(d2));

            // deterministic round-to-nearest: round_ps then truncate
            fs = _mm256_round_ps(fs, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
            __m256i s_int = _mm256_cvttps_epi32(fs);

            // distIC = round(dist*10) + 10000 = round(sqrt(d2)) + 10000
            __m256i dic = _mm256_add_epi32(s_int, i10000);

            // clamp to uint16 range
            dic = _mm256_max_epi32(dic, izero);
            dic = _mm256_min_epi32(dic, i65535);

            // FIXED PACK: pack 8x int32 -> 8x uint16 contiguously
            __m128i dic_lo = _mm256_castsi256_si128(dic);
            __m128i dic_hi = _mm256_extracti128_si256(dic, 1);
            __m128i out8u16 = _mm_packus_epi32(dic_lo, dic_hi);

            _mm_storeu_si128((__m128i*)(outIC + base + t), out8u16);
        }

        // scalar tail
        for (; t < len; ++t) {
            uint32_t j = start + t;

            int32_t dx = (int32_t)xs[j] - (int32_t)xs[i];
            int32_t dy = (int32_t)ys[j] - (int32_t)ys[i];
            int32_t dz = (int32_t)zs[j] - (int32_t)zs[i];
            int32_t d2 = dx*dx + dy*dy + dz*dz;

            float s = std::sqrt((float)d2);                 // deci-units
            uint32_t si = (uint32_t)std::floor(s + 0.5f);   // round-to-nearest for s>=0
            uint32_t dic = si + 10000u;
            if (dic > 65535u) dic = 65535u;

            outIC[base + t] = (uint16_t)dic;
        }
    }
}

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
    const uint L = 10;
    for (uint i = 0; i < L; ++i)
        append(
            float(rand()%100),
            float(rand()%100),
            float(rand()%100));

	vector<uint16_t> ICs;
	ICs.reserve(3*L);
	for (uint32_t i = 0; i < L; ++i)
		{
		ICs.push_back(CoordToIC(x[i]));
		ICs.push_back(CoordToIC(y[i]));
		ICs.push_back(CoordToIC(z[i]));
		}
	BandIndex BI = band_build_index(L, M);
    uint16_t *outIC = myalloc(uint16_t, L*L);
    banded_distances_avx2_u16(ICs.data(), L, M, outIC, BI);
    for (int i = 0; i < L; ++i)
        {
        for (int j = 1; j < i; ++j)
            {
            if (i-j >= M)
                break;
            float d = GetDist(i, j);
            uint16_t dIC = CoordToIC(d);
            uint32_t k = band_ij_to_k(i, j, BI);
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
    //SimpleTest();
    //return;
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
		Chain.GetICs(ICs);

		const uint L = Chain.GetSeqLength();
        const int Li = L;
		BandIndex BI = band_build_index(L, M);
        uint16_t *outIC = myalloc(uint16_t, L*L);
        banded_distances_avx2_u16(ICs.data(), L, M, outIC, BI);

        for (int i = 0; i < Li; ++i)
            {
            for (int j = 1; j < i; ++j)
                {
                if (i-j >= M)
                    break;
                float d = Chain.GetDist(uint(i), uint(j));
                uint16_t dIC = Chain.CoordToIC(d);
                uint k = band_ij_to_k(i, j, BI);
                uint16_t dIC2 = outIC[k];
                fprintf(f, "%d", i);
                fprintf(f, "\t%d", j);
                fprintf(f, "\t%u", dIC);
                fprintf(f, "\t%u", dIC2);
                fprintf(f, "\n");
                }
            }
        break;
		}
	CloseStdioFile(f);
	}
