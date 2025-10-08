#pragma once

#define A(x)	extern atomic<uint> g_ac_##x;
#include "alncounters.h"

static const uint nr_alncounters = 0 
#define A(x)	+ 1
#include "alncounters.h"
;

#define incac(x)	(++g_ac_##x)
#define setac(x, n)	{ g_ac_##x = (n); }

void LogAlnCounts();
