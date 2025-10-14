#include "myutils.h"
#include "features.h"
#include "dss.h"
#include "valuetointtpl.h"

static uint s_DefaultLetters[FEATURE_COUNT];
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
		s_DefaultLetters[i] = 0;
		if (IsFloatFeature((FEATURE) i))
			{
			vector<float> Bins;
			DSS::GetBins((FEATURE) i, s_BinTs[i]);
			}
		}
	return true;
	}
static bool InitDone = Init();

//static uint ValueToInt(float Value, uint AlphaSize, const vector<float> &Ts,
//	uint DefaultLetter)
//	{
//	asserta(DefaultLetter < AlphaSize);
//	if (Value == FLT_MAX)
//		return DefaultLetter;
//
//	asserta(SIZE(Ts) + 1 == AlphaSize);
//	for (uint i = 0; i + 1 < AlphaSize; ++i)
//		if (Value < Ts[i])
//			return i;
//	return AlphaSize - 1;
//	}

uint DSS::ValueToInt_Feature(FEATURE F, float Value)
	{
	assert(uint(F) < FEATURE_COUNT);
	uint AS = DSS::GetAlphaSize(F);
	assert(AS > 0);
	const vector<float> &BinTs = s_BinTs[F];
	uint DefaultLetter = s_DefaultLetters[F];

	if (opt(force_undef))//@@TODO remove this for production
		{
		uint Letter = ValueToIntTpl<true>(Value, AS, BinTs, UINT_MAX);
		assert(Letter < AS || Letter == UINT_MAX);
		return Letter;
		}

// Require 0 <= Letter < AS to allow vector and matrix lookups
//	without special-case testing
	uint Letter = ValueToIntTpl<false>(Value, AS, BinTs, DefaultLetter);
	assert(Letter < AS);
	return Letter;
	}
