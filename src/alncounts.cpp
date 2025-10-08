#include "myutils.h"

#define A(x)	atomic<uint> g_ac_##x;
#include "alncounters.h"

void LogAlnCounts()
	{
	uint SumCounts = 0;
#define A(x)	SumCounts += g_ac_##x;
#include "alncounters.h"
	if (SumCounts == 0)
		return;

#define A(x)	\
	{ uint n = g_ac_##x.load(); \
	Log("%16.16s  %10u", #x, n); \
	if (n > 9999) \
		Log(" (%s)", IntToStr(n)); \
	Log("\n"); \
	}
#include "alncounters.h"
	}
