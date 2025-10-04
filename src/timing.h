#pragma once

#ifndef TIMING
#define TIMING				1
#endif

void InitTiming();
void LogTiming();

#if	TIMING
#include "getticks.h"

enum COUNTER
	{
#define C(x)	COUNTER_##x,
#include "counters.h"
	};

const unsigned CounterCount = 0
#define	C(x) +1
#include "counters.h"
	;

enum TIMER
	{
#define T(x)	TIMER_##x,
#include "timers.h"
	};

const unsigned TimerCount = 0
#define T(x)	+1
#include "timers.h"
	;

const char *TimerToStr(TIMER t);

extern unsigned g_Counters[CounterCount];
#define AddCounter(x, N)	(g_Counters[COUNTER_##x] += (N))
#define IncCounter(x)		(++(g_Counters[COUNTER_##x]))

#define T(x)	extern unsigned g_TimerStartCount##x;
#include "timers.h"

extern TICKS g_LastTicks;
extern TIMER g_CurrTimer;
extern TIMER g_LastTimer;
extern TICKS g_Ticks2D[TimerCount][TimerCount];
extern unsigned g_Counts2D[TimerCount][TimerCount];

void ResetTimers();


#define StartTimer(x)	\
	{	\
	TICKS t = GetClockTicks(); \
	g_CurrTimer = TIMER_##x;	\
	g_Counts2D[g_LastTimer][TIMER_##x] += 1;	\
	g_LastTimer = TIMER_##x;	\
	g_LastTicks = t;	\
	}

#define EndTimer(x)		{ TICKS t = GetClockTicks(); EndTimer_CheckCurr(x); EndTimer_Base(x); g_LastTicks = t; }

struct TimerData
	{
	string Name;
	double Ticks;
	double Calls;
	bool Is2D;
	bool Is1D;

	bool operator <(const TimerData &rhs) const
		{
		return Ticks > rhs.Ticks;
		}
	};

struct GlobalTimingData
	{
	double ElapsedSecs;
	double ElapsedTicks;
	double TimerOverheadTicks;
	};

#else	// TIMING

#define IncCounter(x)		/* empty */
#define AddCounter(x, n)	/* empty */

#define StartTimer(x)		/* empty */
#define PauseTimer(x)		/* empty */
#define ResumeTimer(x)		/* empty */
#define EndTimer(x)			/* empty */

#define StartSTimer(x)		/* empty */
#define PauseSTimer(x)		/* empty */
#define ResumeSTimer(x)		/* empty */
#define EndSTimer(x)		/* empty */

#define ResetTimers()		/* empty */

#endif	// TIMING
