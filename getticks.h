#pragma once

#ifdef _MSC_VER
#include <Windows.h>
typedef unsigned __int64 TICKS;

//#pragma warning(disable:4035)
//inline TICKS GetClockTicks()
//	{
//	_asm
//		{
//		_emit	0x0f
//		_emit	0x31
//		}
//	}

#define	GetClockTicks	__rdtsc

#elif __GNUC__

typedef uint64_t TICKS;
__inline__ uint64_t GetClockTicks()
	{
	uint32_t lo, hi;
	/* We cannot use "=A", since this would use %rax on x86_64 */
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return (uint64_t)hi << 32 | lo;
	}

#else
#error	"getticks_h, unknown compiler"
#endif
