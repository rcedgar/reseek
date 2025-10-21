#include "myutils.h"
#include "features.h"
#include "dss.h"
#include "valuetointtpl.h"

static vector<vector<float> > s_BinTsVec(FEATURE_COUNT);

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
			DSSParams::GetBinTs((FEATURE) i, s_BinTsVec[i]);
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
	const vector<float> &BinTs = s_BinTsVec[F];

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

void DSSParams::OverwriteBinTs(FEATURE F,const vector<float> &BinTs)
	{
	const uint n = SIZE(BinTs);
	asserta(uint(F) < SIZE(s_BinTsVec));
	asserta(SIZE(s_BinTsVec[F]) == n);
	s_BinTsVec[F] = BinTs;
	}
