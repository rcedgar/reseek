#ifndef omplock_h
#define omplock_h

#define TRACE_GLOBAL_LOCKS	0
#define TIME_LOCKS			0

#if	TIME_LOCKS
#include "getticks.h"
void IncLockTicks(TICKS t);
#endif

#if USE_OMP==1

static omp_lock_t g_Lock;

static bool omp_lock_init()
	{
	omp_init_lock(&g_Lock);
	return true;
	}
static bool omp_lock_init_done = omp_lock_init();

static inline void Lock(const char *Msg = "")
	{
#if	TRACE_GLOBAL_LOCKS
	Log("%s:%d %d: Global lock %p %s request\n",
	  __FILE__, __LINE__, omp_get_thread_num(), &g_Lock, Msg);
#endif
#if	TIME_LOCKS
	TICKS t1 = GetClockTicks();
	omp_set_lock(&g_Lock);
	TICKS t2 = GetClockTicks();
	IncLockTicks(t2 - t1);
#else
	omp_set_lock(&g_Lock);
#endif
#if	0 // TRACE_GLOBAL_LOCKS
	Log("%s:%d %d: Global lock %p %s succeed\n",
	  __FILE__, __LINE__, omp_get_thread_num(), &g_Lock, Msg);
#endif
	}

static inline void Unlock(const char *Msg = "")
	{
#if	TRACE_GLOBAL_LOCKS
	Log("%s:%d %d: Global unlock %p %s\n",
	  __FILE__, __LINE__, omp_get_thread_num(), &g_Lock, Msg);
#endif
	omp_unset_lock(&g_Lock);
	}

#endif // USE_OMP

#if USE_OMP==1
#define LOCK(Msg)		Lock(Msg)
#define UNLOCK(Msg)		Unlock(Msg)
#else
#define LOCK(Msg)		/* empty */
#define UNLOCK(Msg)		/* empty */
#endif

#endif // omplock_h
