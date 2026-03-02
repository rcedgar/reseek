#include "myutils.h"
#include "flat_base.h"

atomic<int64_t> g_flat_creates[FE_N];
atomic<int64_t> g_flat_destroys[FE_N];
atomic<int64_t> g_flat_bytes[FE_N];

void log_flat_stats(const string &msg)
	{
	Log("\n");
	Log("log_flat_stats(%s)\n", msg.c_str());
	Log("%10.10s", "Creates");
	Log("  %10.10s", "Destroys");
	Log("  %10.10s", "Actives");
	Log("  %10.10s", "Bytes");
	Log("  Class\n");

#define x(name)	\
	Log("%10lld", g_flat_creates[FE_##name].load()); \
	Log("  %10.10s", Int64ToStr(g_flat_destroys[FE_##name].load())); \
	Log("  %10.10s", Int64ToStr(g_flat_creates[FE_##name].load() - g_flat_destroys[FE_##name].load())); \
	Log("  %10.10s", Int64ToStr(g_flat_bytes[FE_##name].load())); \
	Log("  " #name "\n");
#include "flat_type_names.h"
	}

#if 0
void cmd_test_flat()
	{
	museq_t *museq = create_museq(123);
	log_flat_stats("one");
	down(museq);
	log_flat_stats("two");
	}
#endif // 0