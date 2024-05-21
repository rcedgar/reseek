#pragma once

#include "scop40bench.h"

class Sweeper
	{
public:
	vector<string> m_ParamNames;
	uint m_BestScore = 0;
	uint m_FirstScore = UINT_MAX;
	SCOP40Bench m_SB;
	FILE *m_fFev = 0;

public:
	uint Run(const DSSParams &Params);
	bool Explore1(DSSParams &Params, uint Idx,
	  float Delta, float Z, uint MaxTries);
	};
