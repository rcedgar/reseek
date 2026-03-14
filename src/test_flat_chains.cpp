#include "myutils.h"
#include "dss.h"
#include "chaq.h"
#include "flat_base.h"
#include "flat_chain.h"
#include "pdbchain.h"
#include "flat_distmx.h"
#include "pdbfilescanner.h"
#include "flat_chain_reader.h"

static const uint M = 256;

static bool test_dist_mx(DSS &D, const flat_chain_t* chain)
	{
	const PDBChain &Chain = *D.m_Chain;
	auto dm = chaindistmx_t::newflat(0);
	chaq::create_distmx(chain, dm, M);
	const sid_t *distmx = dm->m_data;
	uint L = D.GetSeqLength();
	asserta(chain->get_length() == L);
	const int Li = L;
	uint counter = 0;
	uint same = 0;
	uint diff1 = 0;
	uint diffgt1 = 0;
	for (int i = 0; i < Li; ++i)
		{
		for (int j = 0; j < i; ++j)
			{
			if (i-j >= int(M))
				continue;
			++counter;
			float d = Chain.GetDist(uint(i), uint(j));
			uint k = banded_ij_to_k(M, i, j);
			asserta(k < L*M);
			sid_t sid = distmx[k];
			float d2 = sid2dist(sid);
			float diff = fabs(d2 - d);
			if (diff < 0.1)
				++same;
			else if (abs(diff) < 1)
				++diff1;
			else
				++diffgt1;
			}
		}
	if (diff1 > 0 || diffgt1 > 0)
		Log("same=%7u (%5.1f%%), diff1=%7u (%5.1f%%), diffgt1=%7u (%5.1f%%) %s\n",
			same, GetPct(same, counter),
			diff1, GetPct(diff1, counter),
			diffgt1, GetPct(diffgt1, counter),
			chain->m_label.c_str());
	chaindistmx_t::release(dm);
	return diffgt1 == 0;
	}

static double test_nn(DSS &D, const flat_chain_t *chain)
	{
	uint L = D.GetSeqLength();
	if (L < 80)
		return 0;
	asserta(chain->get_length() == L);

	chaindistmx_t *dm;
	chaq::create_distmx(chain, dm, M);

	const uint m = 12;
	nnvec_t *nnvec;
	sidvec_t *nndistvec;

	chaq::create_nenvec(dm->m_data, M, L, m, nnvec, nndistvec);
	chaindistmx_t::release(dm);
	auto v = nnvec->m_data;

	uint nsame = 0;
	uint ndiff = 0;
	for (uint i = 0; i < L; ++i)
		{
		uint nn = v[i];
		uint nn2 = UINT_MAX;
		float mindist = FLT_MAX;
		for (uint j = 0; j < L; ++j)
			{
			int diag = abs(int(i) - int(j));
			if (diag < m || diag > int(M))
				continue;
			float d = D.m_Chain->GetDist(i, j);
			if (d < mindist)
				{
				mindist = d;
				nn2 = j;
				}
			}
		if (nn == nn2)
			++nsame;
		else
			++ndiff;
		}
	return float(ndiff)/(nsame + ndiff);
	}

static void test_ic()
	{
	uint same = 0;
	uint diff = 0;
	for (uint32_t ic32 = 0; ic32 < UINT16_MAX; ++ic32)
		{
		ic_t ic = uint16_t(ic32);
		asserta(ic == ic32);
		float d = ic2coord(ic);
		ic_t ic2 = coord2ic(d);
		if (ic == ic2)
			++same;
		else
			++diff;
		}
	if (diff == 0)
		ProgressLog("ic_t <-> float PASS same=%u\n", same);
	else
		ProgressLog("ic_t <-> float DIFFS same=%u, diffs=%u\n", same, diff);

	float resolution = ic2coord(1) - ic2coord(0);
	ProgressLog("Max coord %.1f, resolution %.2f\n",
		ic2coord(UINT16_MAX-1), resolution);
	}

static void test_sid()
	{
	uint same = 0;
	uint diff = 0;
	for (uint32_t sid32 = 0; sid32 < UINT16_MAX; ++sid32)
		{
		sid_t sid = uint16_t(sid32);
		asserta(sid == sid32);
		float d = sid2dist(sid);
		sid_t sid2 = dist2sid(d);
		if (sid == sid2)
			++same;
		else
			++diff;
		}
	if (diff == 0)
		ProgressLog("sid_t <-> float PASS same=%u\n", same);
	else
		Die("sid_t <-> float FAIL same=%u, diffs=%u\n", same, diff);

	float resolution = sid2dist(1) - sid2dist(0);
	ProgressLog("Max sid_t %.1f Angstroms, resolution %.2f A\n",
		sid2dist(UINT16_MAX-1), resolution);
	}

static void test_xyzpair()
	{
	for (uint iter = 0; iter < 1000; ++iter)
		{
		float x1 = float(randu32()%100);
		float y1 = float(randu32()%100);
		float z1 = float(randu32()%100);

		float x2 = float(randu32()%100);
		float y2 = float(randu32()%100);
		float z2 = float(randu32()%100);

		float dx = x1 - x2;
		float dy = y1 - y2;
		float dz = z1 - z2;

		float d = sqrtf(dx*dx + dy*dy + dz*dz);
		sid_t sid1 = dist2sid(d);

		ic_t icx1 = coord2ic(x1);
		ic_t icy1 = coord2ic(y1);
		ic_t icz1 = coord2ic(z1);

		ic_t icx2 = coord2ic(x2);
		ic_t icy2 = coord2ic(y2);
		ic_t icz2 = coord2ic(z2);

		asserta(feq(ic2coord(icx1), x1));
		asserta(feq(ic2coord(icy1), y1));
		asserta(feq(ic2coord(icz1), z1));

		asserta(feq(ic2coord(icx2), x2));
		asserta(feq(ic2coord(icy2), y2));
		asserta(feq(ic2coord(icz2), z2));

		sid_t sid2 = icxyzpair2sid(
			icx1, icy1, icz1,
			icx2, icy2, icz2);

		int diff = abs(int(sid2) - int(sid1));
		double fractdiff = (sid1 == 0 && sid2 == 0 ? 0 : 2.0*double(diff)/(sid1 + sid2));
		asserta(fractdiff < 0.01);
		}
	ProgressLog("test_xyzpair PASS\n");
	}

static void test_neighbor()
	{
/***
Neighbors 3.81 A

>d1v05a_/b.1.18.10
aa         X       Y        Z     icx     icy    icz
A       45.6    26.0    -11.6	10456	10266	9984
M       48.5    24.1    -10.0	10485	10241	9900
***/

	float x1 = 45.6f;
	float x2 = 48.5f;

	float y1 = 26.0f;
	float y2 = 24.1f;

	float z1 = -11.6f;
	float z2 = -10.0f;

	float dx = x1 - x2;
	float dy = y1 - y2;
	float dz = z1 - z2;

	float d = sqrt(dx*dx + dy*dy + dz*dz);
	asserta(feq(d, 3.81));

	ic_t icx1 = coord2ic(x1);
	ic_t icy1 = coord2ic(y1);
	ic_t icz1 = coord2ic(z1);

	ic_t icx2 = coord2ic(x2);
	ic_t icy2 = coord2ic(y2);
	ic_t icz2 = coord2ic(z2);

	// sid=91
	sid_t sid = dist2sid(3.81f);
	asserta(sid == 91);

	// dx = -29, dy = 19, dz = -16
	// (dx*dx + dy*dy + dz*dz) = 1458
	// sid = 1458/16 = 91
	sid_t sid2 = icxyzpair2sid(
		icx1, icy1, icz1,
		icx2, icy2, icz2);
	asserta(sid2 == 91);

	float d2 = sid2dist(sid);
	asserta(feq(d2, 3.81));
	}

static void test_one()
	{
/***
>d1v05a_/b.1.18.10
A	45.6	26.0	-11.6	# 0
M	48.5	24.1	-10.0	# 1
G	49.1	26.3	-7.0	# 2
S	46.8	28.3	-4.8	# 3
D	45.7	31.8	-3.9	# 4
A	43.8	32.3	-0.7	# 5

S	46.8	28.3	-4.8	# 3
A	43.8	32.3	-0.7	# 5

i=3 j=5
 xyz(3) = 10468,10283,9952 = 46.8, 28.3, -4.8
 xyz(5) = 10438,10323,9993 = 43.8, 32.3, -0.7
 = 6.5

sid = 261 = 6.5 A
***/

	float x_i = 46.8f;
	float y_i = 28.3f;
	float z_i = -4.8f;

	float x_j = 43.8f;
	float y_j = 32.3f;
	float z_j = -0.7f;

	float dx = x_i - x_j;
	float dy = y_i - y_j;
	float dz = z_i - z_j;

	float d = sqrtf(dx*dx + dy*dy + dz*dz);
	sid_t sid = dist2sid(d);

	ic_t icx_i = coord2ic(x_i);
	ic_t icy_i = coord2ic(y_i);
	ic_t icz_i = coord2ic(z_i);

	ic_t icx_j = coord2ic(x_j);
	ic_t icy_j = coord2ic(y_j);
	ic_t icz_j = coord2ic(z_j);

	assert(feq(ic2coord(icx_i), x_i));
	assert(feq(ic2coord(icy_i), y_i));
	assert(feq(ic2coord(icz_i), z_i));

	assert(feq(ic2coord(icx_j), x_j));
	assert(feq(ic2coord(icy_j), y_j));
	assert(feq(ic2coord(icz_j), z_j));

	sid_t sid2 = icxyzpair2sid(
		icx_i, icy_i, icz_i,
		icx_j, icy_j, icz_j);
	float d2 = sid2dist(sid2);
	asserta(feq(d, d2));
	}

void cmd_test_flat_dist_types()
	{
	test_one();
	test_neighbor();
	test_xyzpair();
	test_sid();
	test_ic();
	}

void cmd_test_flat_chains()
	{
	vector<PDBChain *> Chains;
	vector<flat_chain_t *> chains;
	ReadChains(g_Arg1, Chains);
	const uint ChainCount = SIZE(Chains);

	DSS D;
	PDBFileScanner FS;
	FS.Open(g_Arg1);
	flat_chain_reader CR;
	CR.Open(FS);
	uint N = 0;
	uint n = 0;
	uint n_fract_diff_gt_1pct = 0;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Processing %.3g%% errs.", GetPct(n, N));
		flat_chain_t* chain = CR.GetNext();
		asserta(chain);
		const PDBChain &Chain = *Chains[ChainIdx];
		asserta(chain->m_label == Chain.m_Label);

		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		if (L < 8)
			continue;

		bool ok = test_dist_mx(D, chain);
		++N;
		if (!ok)
			++n;
		_chkmem();

		double fract_diff = test_nn(D, chain);
		if (fract_diff > 0.01)
			++n_fract_diff_gt_1pct;
		_chkmem();
		}
	ProgressLog("%u / %u (%.2f%%) distmx with diffs > 1\n",
		n, N, GetPct(n, N));
	ProgressLog("%u / %u (%.2f%%) nn with diffs > 1%%\n",
		n_fract_diff_gt_1pct, N, GetPct(n_fract_diff_gt_1pct, N));
	log_flat_stats();
	}
