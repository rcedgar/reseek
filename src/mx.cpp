#include "myutils.h"
#include "mx.h"

mutex g_MxLock;
atomic<int> g_MxNetAllocCount;
atomic<int> g_MxNetAllocCount2;

void LogMxAllocCounts(const char *s)
	{
	ProgressLog("Mx net allocs %s %d %d\n",
				s,
				g_MxNetAllocCount.load(),
				g_MxNetAllocCount2.load());
	}