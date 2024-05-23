#include "myutils.h"

#if USE_OMP
#define L(x)	omp_lock_t g_Lock##x;
#include "lockobjs.h"

static bool Init()
	{
#define L(x)	omp_init_lock(&g_Lock##x);
#include "lockobjs.h"
	return true;
	}
static bool g_InitDone = Init();
#else
#include <mutex>
#define L(x)	mutex g_Lock##x;
#include "lockobjs.h"

#endif
