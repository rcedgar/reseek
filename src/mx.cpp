#include "myutils.h"
#include "mx.h"
#include <mutex>

#define StartTimer(x)	/* empty */
#define EndTimer(x)	/* empty */

uint g_MxAllocCount;
uint g_MxFreeCount;
mutex g_MxLock;

void LogMxAllocCounts(const char *s)
	{
	ProgressLog("Mx allocs %s %u, frees %u, active %u\n",
				s, g_MxAllocCount, g_MxFreeCount, g_MxAllocCount - g_MxFreeCount);
	}