#include "myutils.h"
#include "chaq.h"
#include "quantize.h"

// src/2025-10_reseek_tune/2026-03-25_logodds_and_bins
// [428f078]+ Fix bug in flat_quantize for nen/rendist

static uint16_t ts_mendist16[16-1]={145,178,242,323,403,481,554,631,736,880,1077,1339,1737,2150,3184};
static uint16_t ts_mendist3[3-1]={429,1003};
static uint16_t ts_mendist4[4-1]={323,630,1337};
static uint16_t ts_mendist6[6-1]={215,429,630,1003,1903};
static uint16_t ts_nendist16[16-1]={129,148,168,194,240,289,341,397,458,520,587,677,816,1061,1565};
static uint16_t ts_nendist3[3-1]={256,563};
static uint16_t ts_nendist4[4-1]={194,397,676};
static uint16_t ts_nendist6[6-1]={161,256,397,563,880};
static uint16_t ts_pendist16[16-1]={145,179,247,327,411,490,563,644,757,909,1116,1398,1799,2178,3216};
static uint16_t ts_pendist3[3-1]={439,1041};
static uint16_t ts_pendist4[4-1]={327,645,1399};
static uint16_t ts_pendist6[6-1]={221,439,645,1041,1962};
static uint16_t ts_rendist16[16-1]={308,453,543,617,710,820,949,1101,1276,1489,1748,2026,2199,2594,3634};
static uint16_t ts_rendist3[3-1]={744,1654};
static uint16_t ts_rendist4[4-1]={616,1100,2025};
static uint16_t ts_rendist6[6-1]={516,744,1100,1655,2291};

cp_uint16_t chaq::get_thresholds(FAN fan, uint alpha_size)
	{
#define x(name, size)	if (fan == FAN_##name && alpha_size == size) return ts_##name##size
	x(nendist, 3);
	x(nendist, 4);
	x(nendist, 6);
	x(nendist, 16);

	x(rendist, 3);
	x(rendist, 4);
	x(rendist, 6);
	x(rendist, 16);

	x(pendist, 3);
	x(pendist, 4);
	x(pendist, 6);
	x(pendist, 16);

	x(mendist, 3);
	x(mendist, 4);
	x(mendist, 6);
	x(mendist, 16);
#undef x

	Die("chaq::get_thresholds(%s,%u)", FAN2str(fan), alpha_size);
	return 0;
	}

static uint16_t median_nendist = 397;
static uint16_t median_rendist = 1101;
static uint16_t median_pendist = 645;
static uint16_t median_mendist = 631;

uint16_t chaq::get_undef_value(FAN fan, uint alpha_size)
	{
#define x(name)		if (fan == FAN_##name) return median_##name
	x(nendist);
	x(rendist);
	x(pendist);
	x(mendist);

	Die("get_undef_value(%s)", FAN2str(fan));
	return 0;
	}
