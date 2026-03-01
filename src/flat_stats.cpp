#include "flat_stats.h"
#include <atomic>
#include <cassert>

#ifdef FLAT_TRACK_ALLOCS
    #include <unordered_map>
    #include <mutex>
#endif

namespace flatstats
{
    static std::atomic<uint64_t> g_creates[ft_n];
    static std::atomic<uint64_t> g_deletes[ft_n];
    static std::atomic<uint64_t> g_bytes[ft_n];

#ifdef FLAT_TRACK_ALLOCS
    struct Rec {
        flat_type ft;
        std::size_t bytes;
        const char* type_name;
        const char* file;
        int line;
    };

    static std::unordered_map<void*, Rec> g_live;
    static std::mutex g_live_mx;
#endif

    void on_create(flat_type ft, std::size_t bytes,
                   const char* type_name,
                   const char* file, int line,
                   void* obj_ptr)
    {
        g_creates[ft].fetch_add(1, std::memory_order_relaxed);
        g_bytes[ft].fetch_add((uint64_t)bytes, std::memory_order_relaxed);

#ifdef FLAT_TRACK_ALLOCS
        // file/line may be null if caller didn't pass them; still track.
        std::lock_guard<std::mutex> lock(g_live_mx);
        Rec r;
        r.ft = ft;
        r.bytes = bytes;
        r.type_name = type_name;
        r.file = file;
        r.line = line;
        g_live[obj_ptr] = r;
#endif
    }

    void on_delete(flat_type ft, std::size_t bytes, void* obj_ptr)
    {
        g_deletes[ft].fetch_add(1, std::memory_order_relaxed);
        g_bytes[ft].fetch_sub((uint64_t)bytes, std::memory_order_relaxed);

#ifdef FLAT_TRACK_ALLOCS
        std::lock_guard<std::mutex> lock(g_live_mx);
        std::unordered_map<void*, Rec>::iterator it = g_live.find(obj_ptr);
        if (it != g_live.end())
            g_live.erase(it);
        else
            assert(false && "flatstats: delete of untracked object (double free?)");
#endif
    }

    uint64_t creates(flat_type ft) { return g_creates[ft].load(std::memory_order_relaxed); }
    uint64_t deletes(flat_type ft) { return g_deletes[ft].load(std::memory_order_relaxed); }
    uint64_t bytes(flat_type ft)   { return g_bytes[ft].load(std::memory_order_relaxed); }

    std::size_t report_leaks()
    {
#ifdef FLAT_TRACK_ALLOCS
        std::lock_guard<std::mutex> lock(g_live_mx);

        // Minimal “old style” reporting: you can replace with your logger.
        // Prints one line per live object.
        for (std::unordered_map<void*, Rec>::const_iterator it = g_live.begin();
             it != g_live.end(); ++it)
        {
            const void* p = it->first;
            const Rec& r = it->second;

            // Using stdio to stay “old C++” friendly.
            // (You can swap to iostreams if preferred.)
            std::fprintf(stderr, "LEAK %p ft=%d bytes=%zu type=%s at %s:%d\n",
                         p, (int)r.ft, r.bytes,
                         (r.type_name ? r.type_name : "?"),
                         (r.file ? r.file : "?"), r.line);
        }
        return g_live.size();
#else
        return 0;
#endif
    }
}