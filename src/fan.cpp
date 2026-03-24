#include "myutils.h"
#include "fan.h"

const char *FAN2str(FAN fan)
	{
	switch (fan)
		{
#define f(x)	case FAN_##x: return #x;
#include "flat_alphanamelist.h"
		}
	return "FAN_invalid";
	}

FAN str2FAN(const char *s)
	{
	if (0) ;
#define f(x)	else if (strcmp(s, #x) == 0) return FAN_##x;
#include "flat_alphanamelist.h"
	Die("str2FA(%s)",s);
	return FAN_COUNT;
	}