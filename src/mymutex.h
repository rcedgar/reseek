#pragma once

#include <mutex>
#include "getticks.h"

class mymutex : public std::mutex
	{
public:
	static time_t m_t0;
	static TICKS m_ticks0;

public:
	uint m_nrcalls = 0;
	TICKS m_blocked_ticks = 0;

public:
	void mylock()
		{
		TICKS t1 = GetClockTicks();
		std::mutex::lock();
		m_blocked_ticks += GetClockTicks() - t1;
		++m_nrcalls;
		}

	void myunlock()
		{
		std::mutex::unlock();
		}

	bool mytry_lock()
		{
		++m_nrcalls;
		bool ok = std::mutex::try_lock();
		return ok;
		}

	void logme(const char *name = "mymutex");
	};
