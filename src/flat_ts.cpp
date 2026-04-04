#include "myutils.h"
#include "chaq.h"
#include "quantize.h"

// src/2025-10_reseek_tune/2026-03-25_logodds_and_bins
// [428f078]+ Fix bug in flat_quantize for nen/rendist

// 2025-10_reseek_tune/2026-04-03_logodds_and_bins

static uint16_t ts_mendist3[3-1]={429,1003};
static uint16_t ts_mendist4[4-1]={323,630,1337};
static uint16_t ts_mendist6[6-1]={215,429,630,1003,1903};
static uint16_t ts_mendist8[8-1] = {178,323,481,630,879,1337,2149};
static uint16_t ts_mendist16[16-1]={145,178,242,323,403,481,554,631,736,880,1077,1339,1737,2150,3184};
static uint16_t ts_mendist32[32-1] = {128,145,161,178,203,242,283,323,362,402,442,481,519,554,589,630,679,735,801,878,969,1075,1196,1337,1516,1735,1986,2149,2434,3179,5728};

static uint16_t ts_nendist3[3-1]={256,563};
static uint16_t ts_nendist4[4-1]={194,397,676};
static uint16_t ts_nendist6[6-1]={161,256,397,563,880};
static uint16_t ts_nendist8[8-1] = {148,194,289,397,520,677,1061};
static uint16_t ts_nendist16[16-1]={129,148,168,194,240,289,341,397,458,520,587,677,816,1061,1565};
static uint16_t ts_nendist32[32-1] = {117,129,138,147,157,167,178,192,213,237,262,287,313,339,366,395,425,455,486,518,551,585,625,675,736,814,917,1058,1255,1561,2070};

static uint16_t ts_pendist3[3-1]={439,1041};
static uint16_t ts_pendist4[4-1]={327,645,1399};
static uint16_t ts_pendist6[6-1]={221,439,645,1041,1962};
static uint16_t ts_pendist8[8-1] = {179,327,490,644,909,1398,2178};
static uint16_t ts_pendist16[16-1]={145,179,247,327,411,490,563,644,757,909,1116,1398,1799,2178,3216};
static uint16_t ts_pendist32[32-1] = {128,145,161,179,207,247,286,327,368,410,451,490,527,563,600,644,696,757,827,909,1004,1116,1246,1399,1583,1800,2032,2179,2486,3220,5831};

static uint16_t ts_rendist3[3-1]={744,1654};
static uint16_t ts_rendist4[4-1]={616,1100,2025};
static uint16_t ts_rendist6[6-1]={516,744,1100,1655,2291};
static uint16_t ts_rendist8[8-1] = {453,617,821,1102,1489,2026,2595};
static uint16_t ts_rendist16[16-1]={308,453,543,617,710,820,949,1101,1276,1489,1748,2026,2199,2594,3634};
static uint16_t ts_rendist32[32-1] = {190,308,391,453,502,542,578,616,660,709,762,819,881,948,1021,1100,1184,1275,1374,1487,1611,1746,1891,2025,2112,2198,2352,2593,2978,3633,4935};

static uint16_t ts_pack32[32-1] = {11,14,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,41,43,45,47,50,54};
static uint16_t ts_pack16[16-1] = {14,17,20,22,24,26,28,30,32,34,36,38,41,45,50};
static uint16_t ts_pack8[8-1] = {17,22,26,30,34,39,45};
static uint16_t ts_pack6[6-1] = {19,25,30,35,42};
static uint16_t ts_pack4[4-1] = {22,30,39};
static uint16_t ts_pack3[3-1] = {25,36};

static uint16_t ts_ppack32[32-1] = {1,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,33,36};
static uint16_t ts_ppack16[16-1] = {4,5,7,9,11,12,13,14,16,18,20,22,24,27,31};
static uint16_t ts_ppack8[8-1] = {5,10,13,16,19,22,27};
static uint16_t ts_ppack6[6-1] = {6,12,16,20,25};
static uint16_t ts_ppack4[4-1] = {10,16,22};
static uint16_t ts_ppack3[3-1] = {12,20};

static uint16_t ts_mpack32[32-1] = {1,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,36};
static uint16_t ts_mpack16[16-1] = {4,6,8,10,11,13,15,16,17,18,20,22,24,26,30};
static uint16_t ts_mpack8[8-1] = {5,8,12,17,20,23,27};
static uint16_t ts_mpack6[6-1] = {6,11,16,20,25};
static uint16_t ts_mpack4[4-1] = {8,16,22};
static uint16_t ts_mpack3[3-1] = {11,20};

cp_uint16_t chaq::get_thresholds(FAN fan, uint alpha_size)
	{
#define x(name, size)	if (fan == FAN_##name && alpha_size == size) return ts_##name##size
	x(nendist, 3);
	x(nendist, 4);
	x(nendist, 6);
	x(nendist, 8);
	x(nendist, 16);
	x(nendist, 32);

	x(rendist, 3);
	x(rendist, 4);
	x(rendist, 6);
	x(rendist, 8);
	x(rendist, 16);
	x(rendist, 32);

	x(pendist, 3);
	x(pendist, 4);
	x(pendist, 6);
	x(pendist, 8);
	x(pendist, 16);
	x(pendist, 32);

	x(mendist, 3);
	x(mendist, 4);
	x(mendist, 6);
	x(mendist, 8);
	x(mendist, 16);
	x(mendist, 32);

	x(pack, 3);
	x(pack, 4);
	x(pack, 6);
	x(pack, 8);
	x(pack, 16);
	x(pack, 32);

	x(ppack, 3);
	x(ppack, 4);
	x(ppack, 6);
	x(ppack, 8);
	x(ppack, 16);
	x(ppack, 32);

	x(mpack, 3);
	x(mpack, 4);
	x(mpack, 6);
	x(mpack, 8);
	x(mpack, 16);
	x(mpack, 32);
#undef x

	Die("chaq::get_thresholds(%s,%u)", FAN2str(fan), alpha_size);
	return 0;
	}

// cd $src/2025-10_reseek_tune/2026-04-03_flat_feature_fa
// grep median *.log
static uint16_t median_nendist = 397;
static uint16_t median_rendist = 1101;
static uint16_t median_pendist = 645;
static uint16_t median_mendist = 631;
static uint16_t median_pack = 31;
static uint16_t median_ppack = 16;
static uint16_t median_mpack = 17;

uint16_t chaq::get_undef_value(FAN fan, uint alpha_size)
	{
#define x(name)		if (fan == FAN_##name) return median_##name
	x(nendist);
	x(rendist);
	x(pendist);
	x(mendist);
	x(pack);
	x(ppack);
	x(mpack);

	Die("get_undef_value(%s)", FAN2str(fan));
	return 0;
	}
