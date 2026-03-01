#pragma once

enum FE
	{
#define x(name)	FE_##name,
#include "flat_type_names.h"
	FE_N
	};