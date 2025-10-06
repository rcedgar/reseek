#include "myutils.h"
#include "features.h"
#include "dss.h"

static uint s_DefaultLetters[FEATURE_COUNT];
static UNDEF_BINNING s_UBs[FEATURE_COUNT];
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
		s_UBs[i] = UB_UndefinedIsZeroOverload;
		if (IsFloatFeature((FEATURE) i))
			{
			vector<float> Bins;
			DSS::GetBins((FEATURE) i, s_BinTs[i]);
			}
		}
	return true;
	}
static bool InitDone = Init();

uint DSS::ValueToInt_Feature(FEATURE F, float Value)
	{
	assert(uint(F) < FEATURE_COUNT);
	uint AS = DSS::GetAlphaSize(F);
	assert(AS > 0);
	UNDEF_BINNING UB = s_UBs[F];
	const vector<float> &BinTs = s_BinTs[F];
	uint DefaultLetter = s_DefaultLetters[F];
	uint Letter = ValueToInt(Value, UB, AS, BinTs, DefaultLetter);
	return Letter;
	}
