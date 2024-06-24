#pragma once

enum USERFIELD
	{
	UF_Undefined,
#define x(name)	UF_##name,
#include "userfieldnames.h"
	};

USERFIELD StrToUF(const char *Str);
USERFIELD StrToUF(const string &Str);
const char *UFToStr(USERFIELD UF);
