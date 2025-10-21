#include "myutils.h"
#include "features.h"
#include "dss.h"
#include "valuetointtpl.h"

static vector<vector<float> > s_BinTs(FEATURE_COUNT);

static bool IsFloatFeature(FEATURE Feat)
	{
	switch (Feat)
		{
#define F(x)	case FEATURE_##x:	return true;
#include "floatfeatures.h"
		}
	return false;
	}

static bool Init()
	{
	for (uint i = 0; i < FEATURE_COUNT; ++i)
		{
		if (IsFloatFeature((FEATURE) i))
			{
			vector<float> Bins;
			DSSParams::GetBins((FEATURE) i, s_BinTs[i]);
			}
		}
	return true;
	}
static bool InitDone = Init();

uint DSSParams::ValueToInt_Feature(FEATURE F, float Value)
	{
	assert(uint(F) < FEATURE_COUNT);
	uint AS = DSSParams::GetAlphaSize(F);
	assert(AS > 0);
	const vector<float> &BinTs = s_BinTs[F];

	if (opt(force_undef))//@@TODO remove this for production
		{
		uint Letter = ValueToIntTpl<true>(Value, AS, BinTs, UINT_MAX);
		assert(Letter < AS || Letter == UINT_MAX);
		return Letter;
		}

// Require 0 <= Letter < AS to allow vector and matrix lookups
//	without special-case testing
	uint Letter = ValueToIntTpl<false>(Value, AS, BinTs, 0);
	assert(Letter < AS);
	return Letter;
	}
