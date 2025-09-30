#if 0 // @@ DELETE
#include "myutils.h"
#include "mx.h"

#if 0
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
#endif // 0
#endif