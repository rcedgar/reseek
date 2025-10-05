#pragma once

#ifndef TIMING
#define TIMING				0
#endif

void LogTiming();

#if	TIMING
#include "getticks.h"

enum COUNTER
	{
#define C(x)	COUNTER_##x,
#include "counters.h"
	};

const uint CounterCount = 0
#define	C(x) +1
#include "counters.h"
	;

enum TIMER
	{
#define T(x)	TIMER_##x,
#include "timers.h"
	};

const uint TimerCount = 0
#define T(x)	+1
#include "timers.h"
	;

const char *TimerToStr(TIMER t);

extern uint g_Counters[CounterCount];
#define AddCounter(x, N)	(g_Counters[COUNTER_##x] += (N))
#define IncCounter(x)		(++(g_Counters[COUNTER_##x]))

#define T(x)	extern uint g_TimerStartCount##x;
#include "timers.h"

extern TICKS g_LastTicks;
extern TIMER g_CurrTimer;
extern TIMER g_LastTimer;
extern TICKS g_Ticks1D[TimerCount];
extern uint g_Counts1D[TimerCount];
extern TICKS g_Ticks2D[TimerCount][TimerCount];
extern uint g_Counts2D[TimerCount][TimerCount];

void ResetTimers();

#define StartTimer(x)	\
	{	\
	if (g_CurrTimer != TIMER_None)	\
		{	\
		const char *s = TimerToStr(g_CurrTimer);	\
		Die("StartTimer(" #x "), curr=%s %d %d", s, g_CurrTimer, TIMER_##x);	\
		}	\
	TICKS t = GetClockTicks(); \
	g_Ticks2D[g_LastTimer][TIMER_##x] += t - g_LastTicks;	\
	g_Counts2D[g_LastTimer][TIMER_##x] += 1;	\
	g_CurrTimer = TIMER_##x;	\
	g_LastTimer = TIMER_##x;	\
	g_LastTicks = t;	\
	}

#define EndTimer(x)	\
	{	\
	if (g_CurrTimer != TIMER_##x)	\
		Die("EndTimer(" #x "), curr=%s %d %d", TimerToStr(g_CurrTimer), TIMER_##x, g_CurrTimer);	\
	TICKS t = GetClockTicks(); \
	g_Ticks1D[TIMER_##x] += t - g_LastTicks;	\
	g_Counts1D[TIMER_##x] += 1;	\
	g_LastTimer = TIMER_##x;	\
	g_LastTicks = t;	\
	g_CurrTimer = TIMER_None;	\
	}

struct TimerData
	{
	int Timer1;
	int Timer2;
	double Ticks;
	double Calls;
	bool Is2D;

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

#define ResetTimers()		/* empty */

#endif	// TIMING
