#include "myutils.h"

#define TEST	0

static uint g_Primes[] =
	{
#include "primes.h"
	};
static uint g_PrimeCount = uint(sizeof(g_Primes)/sizeof(g_Primes[0]));

uint FindPrime(uint Min, uint Max)
	{
	for (unsigned i = 0; i < g_PrimeCount; ++i)
		{
		uint Prime = g_Primes[i];
		if (Prime >= Min && Prime <= Max)
			return Prime;
		}
	return (Min + Max)/2 + 1;
	}
