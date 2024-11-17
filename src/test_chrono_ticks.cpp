#include "myutils.h"
#include "getticks.h"
#include <chrono>

static double Slow(uint32_t n)
    {
    double Sum = 0;
    for (uint i = 0; i < n*100000; ++i)
        Sum += sin(randu32()%100);
    return Sum;
    }

/***
ticks 10688500, cticks 34066829, r=0.314 | MSVC
ticks 12162862, cticks 38765634, r=0.314 | Linux gcc x86
***/
void cmd_test_chrono_ticks()
	{
    double Sum = 0;
    for (uint Try = 0; Try < 10; ++Try)
        {
        uint64_t t0 = GetClockTicks();
        auto c0 = std::chrono::high_resolution_clock::now();

        Sum += Slow(randu32()%10);

        auto c1 = std::chrono::high_resolution_clock::now();
        uint64_t t1 = GetClockTicks();

        uint32_t ticks = uint32_t(t1 - t0);

        auto elapsed = c1 - c0;
        uint32_t cticks = (uint32_t) elapsed.count();
        double r = double(cticks)/ticks;
        ProgressLog("ticks %u, cticks %u, r=%.3g\n",
          cticks, ticks, r);
        }

    Log("Sum=%.3g\n", Sum);
	}
