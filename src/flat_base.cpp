#include "myutils.h"
#include "flat_base.h"

#if TRACK_ACTIVE
atomic<int64_t> g_flat_creates[FE_N];
atomic<int64_t> g_flat_destroys[FE_N];
atomic<int64_t> g_flat_bytes[FE_N];
#endif

#if TRACK_SRC
list<void *> g_flat_obj_list;
mutex g_flat_obj_list_lock;
#endif

void log_flat_stats(const string &msg)
	{
#if TRACK_ACTIVE
	Log("\n");
	Log("log_flat_stats(%s)\n", msg.c_str());
	Log("%10.10s", "Creates");
	Log("  %10.10s", "Destroys");
	Log("  %10.10s", "Actives");
	Log("  %10.10s", "Bytes");
	Log("  Class\n");

	uint64_t total_active = 0;
#define x(name)	\
	{ \
	auto n = g_flat_creates[FE_##name].load(); \
	total_active += g_flat_creates[FE_##name].load() - g_flat_destroys[FE_##name].load(); \
	if (n > 0) \
	{ \
	Log("%10.10s", Int64ToStr(g_flat_creates[FE_##name].load())); \
	Log("  %10.10s", Int64ToStr(g_flat_destroys[FE_##name].load())); \
	Log("  %10.10s", Int64ToStr(g_flat_creates[FE_##name].load() - g_flat_destroys[FE_##name].load())); \
	Log("  %10.10s", Int64ToStr(g_flat_bytes[FE_##name].load())); \
	Log("  " #name "\n"); \
	} }
#include "flat_type_names.h"

	Log("%10.10s  %10.10s  %10.10s  %10.10s  TOTAL ACTIVE\n",
		"", "", Int64ToStr(total_active), "");
#endif // TRACK_ACTIVE

#if TRACK_SRC
	vector<string> srcfiles;
	vector<int> srclines;
	vector<FE> flat_types;
	vector<uint32_t> sizes;
	map<pair<string, int>, uint64_t> loc2totalsize;
	for (list<void *>::const_iterator iter = g_flat_obj_list.begin();
		iter != g_flat_obj_list.end(); ++iter)
		{
		flat_base<uint32_t, FE_BASE> *pobj =
			(flat_base<uint32_t, FE_BASE> *) *iter;
		FE fe = pobj->m_fe;
		string srcfile = "src???";
		int srcline = 0;
		if (pobj->m_srcfile)
			{
			srcfile = string(pobj->m_srcfile);
			srcline = pobj->m_srcline;
			}
		vector<string> flds;
		Split(srcfile, flds, '\\');
		srcfile = flds[flds.size()-1];
		srcfiles.push_back(srcfile);
		srclines.push_back(srcline);
		flat_types.push_back(pobj->m_fe);
		sizes.push_back(pobj->m_size);

		pair<string, int> loc(srcfile, srcline);
		if (loc2totalsize.find(loc) == loc2totalsize.end())
			loc2totalsize[loc] = 0;
		loc2totalsize[loc] += pobj->m_size;

		//Log("%s [%u] %s(%d)\n", FE2str(pobj->m_fe), pobj->m_size,
		//	srcfile.c_str(), srcline);
		}
	std::vector<std::pair<std::pair<std::string,int>, uint64_t>> v(
		loc2totalsize.begin(), loc2totalsize.end());

	std::sort(v.begin(), v.end(),
		[](const auto &a, const auto &b)
		{
			return a.second > b.second;   // decreasing size
		});

	Log("\n");
	Log("%10.10s", "Bytes");
	Log("  Source");
	Log("\n");
	uint64_t total_size = 0;
	for (const auto &e : v)
		{
		const std::string &file = e.first.first;
		int line = e.first.second;
		uint64_t size = e.second;
		total_size += size;

		Log("%10.10s", MemBytesToStr(double(size)));
		Log("  %s(%d)", file.c_str(), line);
		Log("\n");
		}
	Log("%10.10s  TOTAL\n", Int64ToStr(total_size));
#endif
	}

#if 0
void cmd_test_flat()
	{
	auto p1 = chainaa_t::newflat(123);
	auto p2 = chaindistmx_t::newflat(456, 64);

	log_flat_stats();
	}
#endif // 0