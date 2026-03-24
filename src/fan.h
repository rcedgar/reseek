#pragma once

enum FAN
	{
#define f(x)	FAN_##x,
#include "flat_alphanamelist.h"
	FAN_COUNT
	};

const char *FAN2str(FAN fan);
FAN str2FAN(const char *s);