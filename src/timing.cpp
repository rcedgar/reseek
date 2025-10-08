#include "myutils.h"
#include "timing.h"
#include "getticks.h"
#include <time.h>
#include <algorithm>

#ifdef _MSC_VER
#include <Windows.h>	// for IsDebuggerPresent()
#endif

#if	TIMING

TIMER g_CurrTimer = TIMER_None;
TIMER g_LastTimer = TIMER_None;
TICKS g_LastTicks = GetClockTicks();
TICKS g_StartTicks = GetClockTicks();
time_t g_StartTime = time(0);

TICKS g_Ticks1D[TimerCount];
uint g_Counts1D[TimerCount];
TICKS g_Ticks2D[TimerCount][TimerCount];
uint g_Counts2D[TimerCount][TimerCount];

void ResetTimers()
	{
	for (uint i = 0; i < TimerCount; ++i)
		{
		for (uint j = 0; j < TimerCount; ++j)
			{
			g_Ticks2D[i][j] = 0;
			g_Counts2D[i][j] = 0;
			}
		}
	}

static double GetTimerOverheadTicks()
	{
	const uint N = 1000*1000;

	TICKS t1 = GetClockTicks();
	for (uint i = 0; i < N; ++i)
		{
		StartTimer(GetTimerOverhead);
		EndTimer(GetTimerOverhead);
		}
	TICKS t2 = GetClockTicks();

	double Ticks = double(t2 - t1);
	return Ticks/N;
	}

const char *TimerToStr(TIMER t)
	{
	switch (t)
		{
#define T(x)	case TIMER_##x: return #x;
#include "timers.h"
#undef T
		}
	return "?";
	}

static void LogTimers(const GlobalTimingData &GTD,
  const vector<TimerData> &TDs)
	{
	const uint TDCount = SIZE(TDs);

	double TotalTimerCalls = 0.0;
	double TotalTimerTicks = 0.0;
	for (uint i = 0; i < TDCount; ++i)
		{
		const TimerData &TD = TDs[i];
		TotalTimerTicks += TD.Ticks;
		}

	double TicksPerSec = GTD.ElapsedTicks/GTD.ElapsedSecs;
	double TotalTimerSecs = TotalTimerTicks/TicksPerSec;
	double TotalTimerOverheadTicks = TotalTimerCalls*GTD.TimerOverheadTicks;
	double TotalTimerOverheadSecs = TotalTimerOverheadTicks/TicksPerSec;

	double TotalPct = 0.0;
	double TotalPctShown = 0.0;
	double TotalCPct = 0.0;
	double TotalSecs = 0.0;
	double MinPct = opt(min_timer_pct);

	Log("\n");
	Log("MinPct %.1f%%\n", MinPct);
	Log("    Pct   TotPct        Ticks        Secs       Calls  Ticks/Call  Timer\n");
	for (uint i = 0; i < TDCount; ++i)
		{
		const TimerData &TD = TDs[i];

		double pct = GetPct(TD.Ticks, TotalTimerTicks);
		TotalPct += pct;

		double Secs = TD.Ticks/TicksPerSec;
		TotalSecs += Secs;

		double CTicks = TD.Ticks;
		CTicks -= TD.Calls*GTD.TimerOverheadTicks;
		double TicksPerCall = TD.Calls == 0.0 ? 0.0 : CTicks/TD.Calls;

		if (pct >= MinPct)
			{
			TotalPctShown += pct;

			string Name;
			if (TD.Is2D)
				{
				Name = TimerToStr((TIMER) TD.Timer1);
				Name += " - ";
				Name += TimerToStr((TIMER) TD.Timer2);
				}
			else
				Name = TimerToStr((TIMER) TD.Timer1);

			Log("%6.1f%%", pct);
			Log("  %6.1f%%", TotalPct);
			Log("  %11.4e", TD.Ticks);
			Log("  %10.10s", FloatToStr(Secs));
			Log("  %10.10s", FloatToStr(TD.Calls));
			if (TD.Is2D)
				Log("  %10.10s", "");
			else
				Log("  %10.10s", FloatToStr(TicksPerCall));
			Log("  %s", Name.c_str());
			Log("\n");
			}
		}

	Log("=======  =======  ===========  ==========  ==========  ==========  =====\n");

	Log("%6.1f%%", TotalPctShown);
	Log("  %6.1f%%", TotalPct);
	Log("  %11.11s", FloatToStr(TotalTimerTicks));
	Log("  %10.10s", FloatToStr(TotalSecs));
	Log("  %10.10s", FloatToStr(TotalTimerCalls));
	Log("  %10.10s", "");
	Log("  TOTAL\n");
	Log("\n");

	Log("\n");
	Log("%11.4g  Ticks elapsed (%.1f secs)\n",
	  GTD.ElapsedTicks, GTD.ElapsedTicks/TicksPerSec);
	Log("%11.4g  Ticks total timers (%.1f%%)\n",
	  TotalTimerTicks,
	  GetPct(TotalTimerTicks, GTD.ElapsedTicks));
	Log("%11.11s  Secs elapsed\n", SecsToStr(GTD.ElapsedSecs));
	Log("%11.11s  Timer overhead ticks\n", FloatToStr(GTD.TimerOverheadTicks));

	double TimerOverheadSecs = TotalTimerCalls*GTD.TimerOverheadTicks/TicksPerSec;
	double EstdSecs = GTD.ElapsedSecs - TimerOverheadSecs;
	Log("%11.11s  Secs timer overhead (%.1f%%)\n",
	  SecsToStr(TimerOverheadSecs),
	  GetPct(TimerOverheadSecs, GTD.ElapsedSecs));
	Log("%11.11s  Secs estd. without timers\n",
	  SecsToStr(EstdSecs));
	}

void LogTiming()
	{
	TICKS ExitTicks  = GetClockTicks();
	time_t ExitTime = time(0);

	g_Ticks2D[g_LastTimer][TIMER_ExitTiming] = ExitTicks - g_LastTicks;

	GlobalTimingData GTD;
	double ElapsedSecs = double(ExitTime - g_StartTime);

	if (ElapsedSecs == 0.0)
		ElapsedSecs = 1.0;

	double ElapsedTicks = double(ExitTicks - g_StartTicks);

	vector<TimerData> TDs;
	for (uint i = 0; i < TimerCount; ++i)
		{
		double Ticks = (double) g_Ticks1D[i];
		if (Ticks > 0.0)
			{
			TimerData TD;
			TD.Timer1 = i;
			TD.Timer2 = i;
			TD.Is2D = false;
			TD.Ticks = double(g_Ticks1D[i]);
			TD.Calls = g_Counts1D[i];
			TDs.push_back(TD);
			}
		}

	for (uint i = 0; i < TimerCount; ++i)
		for (uint j = 0; j < TimerCount; ++j)
			{
			double Ticks = (double) g_Ticks2D[i][j];
			if (Ticks > 0.0)
				{
				TimerData TD;
				TD.Timer1 = i;
				TD.Timer2 = j;
				TD.Is2D = true;
				TD.Ticks = double(g_Ticks2D[i][j]);
				TD.Calls = g_Counts2D[i][j];
				TDs.push_back(TD);
				}
			}
	sort(TDs.begin(), TDs.end());

	Log("\n");
#ifdef _MSC_VER
	Log("Debugger %s\n", IsDebuggerPresent() ? "Yes " : "No");
#endif

	GTD.ElapsedSecs = ElapsedSecs;
	GTD.ElapsedTicks = ElapsedTicks;
	GTD.TimerOverheadTicks = GetTimerOverheadTicks();

	LogTimers(GTD, TDs);
	}

#else // #if TIMING

void InitTiming()
	{
	}

void LogTiming()
	{
	}

#endif // #if TIMING
