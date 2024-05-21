#include "myutils.h"
#include "features.h"

const char *FeatureToStr(uint FeatureIndex)
	{
	return FeatureToStr(FEATURE(FeatureIndex));
	}

const char *FeatureToStr(FEATURE f)
	{
	switch (f)
		{
#define F(x)	case FEATURE_##x: return #x;
#include "featurelist.h"
		}
	Die("FeatureToStr(%u)", f);
	return "?";
	}

FEATURE StrToFeature(const char *s)
	{
#define F(x)	if (!stricmp(s, #x)) return FEATURE_##x;
#include "featurelist.h"
	Die("StrToFeature(%s)", s);
	return FEATURE(UINT_MAX);
	}

uint StrToFeatureIndex(const char *s)
	{
#define F(x)	if (!stricmp(s, #x)) return uint(FEATURE_##x);
#include "featurelist.h"
	Die("StrToFeatureIndex(%s)", s);
	return UINT_MAX;
	}

bool FeatureIsInt(FEATURE f)
	{
	switch (f)
		{
	case FEATURE_AA: return true;

#define F(x)	case FEATURE_##x: return true;
#include "intfeatures.h"
#undef F

#define F(x)	case FEATURE_##x: return false;
#include "floatfeatures.h"
#undef F
		}
	asserta(false);
	return false;
	}

bool FeatureIsInt(uint FeatureIndex)
	{
	return FeatureIsInt(FEATURE(FeatureIndex));
	}

void cmd_feature_stats()
	{
	for (uint F = 0; F < FEATURE_COUNT; ++F)
		{
		uint AS = g_AlphaSizes2[F];
		float **Mx = g_ScoreMxs2[F];
		ProgressLog("[%2u]  %s",
		  F, FeatureToStr(F));
		if (Mx == 0)
			ProgressLog("  < missing scoremx");
		ProgressLog("\n");
		}
	}

extern float **g_ScoreMxs2[FEATURE_COUNT];
extern uint g_AlphaSizes2[FEATURE_COUNT];
