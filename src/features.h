#pragma once

enum FEATURE
	{
#define F(x)	FEATURE_##x,
#include "featurelist.h"
	};

const uint FEATURE_COUNT = 0 +
#define F(x)	+ 1
#include "featurelist.h"
;

const char *FeatureToStr(uint FeatureIndex);
const char *FeatureToStr(FEATURE f);
FEATURE StrToFeature(const char *s);
uint StrToFeatureIndex(const char *s);
bool FeatureIsInt(FEATURE f);
bool FeatureIsInt(uint FeatureIndex);

//extern float **g_ScoreMxs2[FEATURE_COUNT];
//extern uint g_AlphaSizes2[FEATURE_COUNT];
