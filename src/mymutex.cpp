#include "myutils.h"
#include "mymutex.h"

time_t mymutex::m_t0 = time(0);
TICKS mymutex::m_ticks0 = GetClockTicks();

void mymutex::logme(const char *name)
	{
	time_t t1 = time(0);
	TICKS ticks1 = GetClockTicks();

	time_t secs = t1 - m_t0;
	if (secs == 0)
		secs = 1;
	double ticks_per_sec = double(ticks1)/secs;

	ProgressLog("%s: %u calls, blocked %.3g ticks, %.2f secs (%.3g ticks/sec)\n",
		name,
		m_nrcalls,
		(double) m_blocked_ticks,
		m_blocked_ticks/ticks_per_sec,
		ticks_per_sec);
	}
