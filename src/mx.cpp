#include "myutils.h"
#include "mx.h"

uint g_MxAllocCount;
uint g_MxFreeCount;
mutex g_MxLock;

void LogMxAllocCounts(const char *s)
	{
	ProgressLog("Mx allocs %s %u, frees %u, active %u\n",
				s, g_MxAllocCount, g_MxFreeCount, g_MxAllocCount - g_MxFreeCount);
	}