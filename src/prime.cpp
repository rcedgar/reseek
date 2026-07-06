#include "myutils.h"

#define TEST	0

static uint g_Primes[] =
	{
#include "primes.h"
	};
static uint g_PrimeCount = uint(sizeof(g_Primes)/sizeof(g_Primes[0]));

uint FindPrime(uint Min, uint Max)
	{
	if (Min > Max)
		return (Min + Max)/2 + 1;
	uint lo = 0;
	uint hi = g_PrimeCount;
	while (lo < hi)
		{
		uint mid = (lo + hi)/2;
		if (g_Primes[mid] < Min)
			lo = mid + 1;
		else
			hi = mid;
		}
	if (lo < g_PrimeCount && g_Primes[lo] <= Max)
		return g_Primes[lo];
	return (Min + Max)/2 + 1;
	}
