#pragma once
#include <stdio.h>

#ifdef _WIN32

#include <windows.h>
#include <psapi.h>

struct RunStats
{
    LARGE_INTEGER t0{}, freq{};
};

static double filetime_sec(const FILETIME& ft)
{
    ULARGE_INTEGER u;
    u.LowPart = ft.dwLowDateTime;
    u.HighPart = ft.dwHighDateTime;
    return double(u.QuadPart) * 1e-7;
}

static DWORD get_physical_cores()
{
    DWORD len = 0;
    GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &len);
    if (len == 0) return 0;

    char* buf = (char*)malloc(len);
    if (!buf) return 0;

    if (!GetLogicalProcessorInformationEx(RelationProcessorCore,
                                          (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX)buf,
                                          &len))
    {
        free(buf);
        return 0;
    }

    DWORD cores = 0;
    char* p = buf;
    while (p < buf + len)
    {
        auto info = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX)p;
        if (info->Relationship == RelationProcessorCore)
            ++cores;
        p += info->Size;
    }

    free(buf);
    return cores;
}

static void start_runstats(RunStats& rs)
{
    QueryPerformanceFrequency(&rs.freq);
    QueryPerformanceCounter(&rs.t0);
}

static void log_runstats(const RunStats& rs, FILE* f)
{
    if (!f) return;

    PROCESS_MEMORY_COUNTERS pmc{};
    pmc.cb = sizeof(pmc);
    if (!GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return;

    FILETIME create_ft{}, exit_ft{}, kernel_ft{}, user_ft{};
    if (!GetProcessTimes(GetCurrentProcess(),
        &create_ft, &exit_ft, &kernel_ft, &user_ft))
        return;

    LARGE_INTEGER t1;
    QueryPerformanceCounter(&t1);

    double wall_sec =
        double(t1.QuadPart - rs.t0.QuadPart) / double(rs.freq.QuadPart);

    double user_sec = filetime_sec(user_ft);
    double sys_sec  = filetime_sec(kernel_ft);
    double cpu_sec  = user_sec + sys_sec;

    SYSTEM_INFO si{};
    GetSystemInfo(&si);
    DWORD logical = si.dwNumberOfProcessors;
    if (logical == 0) logical = 1;

    DWORD physical = get_physical_cores();

    double denom = double(logical);

    // Apply SMT correction only if clearly 2x
    if (physical > 0 && logical >= 2 * physical - 1 && logical <= 2 * physical + 1)
        denom = double(physical);

    double cpu_frac = 0.0;
    if (wall_sec > 0)
        cpu_frac = cpu_sec / (wall_sec * denom);

    double peak_mem_bytes = double(pmc.PeakWorkingSetSize);

    //fprintf(f,
    //    "peak_rss_mb=%.1f cpu_used_pct=%.1f user=%.3f sys=%.3f wall=%.3f ncpu_eff=%.0f\n",
    //    double(pmc.PeakWorkingSetSize) / (1024.0 * 1024.0),
    //    100.0 * cpu_frac,
    //    user_sec,
    //    sys_sec,
    //    wall_sec,
    //    denom);

    fprintf(f, "Peak mem %s", MemBytesToStr(peak_mem_bytes));
    fprintf(f, ", CPU %.1f%%", 100.0*cpu_frac);
    fprintf(f, "\n");
}

#else // Linux

#include <sys/resource.h>
#include <unistd.h>
#include <time.h>
#include <set>
#include <utility>
#include <fstream>
#include <string>

struct RunStats
{
    double wall_start = 0.0;
};

static double now_sec()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return double(ts.tv_sec) + 1e-9 * double(ts.tv_nsec);
}

static long get_physical_cores()
{
    std::ifstream in("/proc/cpuinfo");
    if (!in) return 0;

    std::set<std::pair<int,int>> cores;
    std::string line;
    int phys = -1, core = -1;

    while (std::getline(in, line))
    {
        if (line.find("physical id") != std::string::npos)
            phys = atoi(line.c_str() + line.find(":") + 1);
        else if (line.find("core id") != std::string::npos)
            core = atoi(line.c_str() + line.find(":") + 1);
        else if (line.empty())
        {
            if (phys >= 0 && core >= 0)
                cores.insert({phys, core});
            phys = core = -1;
        }
    }

    return cores.empty() ? 0 : (long)cores.size();
}

static void start_runstats(RunStats& rs)
{
    rs.wall_start = now_sec();
}

static void log_runstats(const RunStats& rs, FILE* f)
{
    if (!f) return;

    struct rusage ru{};
    if (getrusage(RUSAGE_SELF, &ru) != 0)
        return;

    double user_sec =
        double(ru.ru_utime.tv_sec) + 1e-6 * double(ru.ru_utime.tv_usec);
    double sys_sec =
        double(ru.ru_stime.tv_sec) + 1e-6 * double(ru.ru_stime.tv_usec);
    double cpu_sec = user_sec + sys_sec;

    double wall_sec = now_sec() - rs.wall_start;

    long logical = sysconf(_SC_NPROCESSORS_ONLN);
    if (logical <= 0) logical = 1;

    long physical = get_physical_cores();

    double denom = double(logical);

    // Apply SMT correction only if ~2x
    if (physical > 0 && logical >= 2 * physical - 1 && logical <= 2 * physical + 1)
        denom = double(physical);

    double cpu_frac = 0.0;
    if (wall_sec > 0)
        cpu_frac = cpu_sec / (wall_sec * denom);

    //fprintf(f,
    //    "peak_rss_mb=%.1f cpu_used_pct=%.1f user=%.3f sys=%.3f wall=%.3f ncpu_eff=%.0f\n",
    //    double(ru.ru_maxrss) / 1024.0,
    //    100.0 * cpu_frac,
    //    user_sec,
    //    sys_sec,
    //    wall_sec,
    //    denom);
    fprintf(f, "Peak mem %s", MemBytesToStr(ru.ru_maxrss*1024.0));
    fprintf(f, ", CPU %.1f%%", 100.0*cpu_frac);
    fprintf(f, "\n");
}

#endif