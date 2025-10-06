#include "myutils.h"
#include "undef_binning.h"
#include "dss.h"

/***
Undefined values:
	Letter is zero.
	DSS:GetFloat_F returns FLT_MAX.

Alphabet has AlphaSize (AS) letters.
Letter 0 is undefined.
That leaves AS-1 defined letters separated by AS-2 thresholds.
***/

uint DSS::GetBinThresholdCount(uint AlphaSize, UNDEF_BINNING UB)
	{
	switch (UB)
		{
	case UB_NeverUndefined:				return AlphaSize - 1;
	case UB_UndefinedIsDefaultLetter:	return AlphaSize - 1;
	case UB_IgnoreUndefined:			return AlphaSize - 1;
	case UB_UndefinedIsZeroOverload:	return AlphaSize - 1;
	case UB_UndefinedIsOnlyZero:		return AlphaSize - 2;
		}
	Die("GetBinThresholdCount(%u, %s)", AlphaSize, UBToStr(UB));
	return UINT_MAX;
	}

uint DSS::ValueToInt(float Value, UNDEF_BINNING UB, uint AlphaSize,
					 const vector<float> &Ts, uint DefaultLetter)
	{
	uint Letter = UINT_MAX;
	switch (UB)
		{
	case UB_NeverUndefined:
		Letter = ValueToInt_Never(Value, AlphaSize, Ts, DefaultLetter);
		break;

	case UB_UndefinedIsOnlyZero:
		Letter = ValueToInt_OnlyZero(Value, AlphaSize, Ts, DefaultLetter);
		break;

	case UB_UndefinedIsZeroOverload:
		Letter = ValueToInt_ZeroOverload(Value, AlphaSize, Ts, DefaultLetter);
		break;

	case UB_UndefinedIsDefaultLetter:
		Letter = ValueToInt_Default(Value, AlphaSize, Ts, DefaultLetter);
		break;

	case UB_IgnoreUndefined:
		Letter = ValueToInt_Ignore(Value, AlphaSize, Ts, DefaultLetter);
		break;

	default:
		asserta(false);
		}

	asserta(Letter < AlphaSize);
	return Letter;
	}

uint DSS::ValueToInt_Never(float Value, uint AlphaSize,
						   const vector<float> &Ts, uint DefaultLetter)
	{
	asserta(Value != FLT_MAX);
	asserta(DefaultLetter == UINT_MAX);
	asserta(SIZE(Ts) + 1 == AlphaSize);
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		if (Value < Ts[i])
			return i;
	return AlphaSize - 1;
	}

uint DSS::ValueToInt_Ignore(float Value, uint AlphaSize,
						   const vector<float> &Ts, uint DefaultLetter)
	{
	asserta(Value != FLT_MAX);
	asserta(DefaultLetter == UINT_MAX);
	asserta(SIZE(Ts) + 1 == AlphaSize);
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		if (Value < Ts[i])
			return i;
	return AlphaSize - 1;
	}

uint DSS::ValueToInt_OnlyZero(float Value, uint AlphaSize,
						   const vector<float> &Ts, uint DefaultLetter)
	{
	asserta(DefaultLetter == 0);
	if (Value == FLT_MAX)
		return 0;

	asserta(SIZE(Ts) + 2 == AlphaSize);
	for (uint i = 0; i + 2 < AlphaSize; ++i)
		if (Value < Ts[i])
			return i + 1;
	return AlphaSize - 1;
	}

uint DSS::ValueToInt_Default(float Value, uint AlphaSize,
						   const vector<float> &Ts, uint DefaultLetter)
	{
	asserta(DefaultLetter < AlphaSize);
	if (Value == FLT_MAX)
		return DefaultLetter;

	asserta(SIZE(Ts) + 1 == AlphaSize);
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		if (Value < Ts[i])
			return i;
	return AlphaSize - 1;
	}

uint DSS::ValueToInt_ZeroOverload(float Value, uint AlphaSize,
						   const vector<float> &Ts, uint DefaultLetter)
	{
	asserta(DefaultLetter == 0);
	if (Value == FLT_MAX)
		return 0;

	asserta(SIZE(Ts) + 1 == AlphaSize);
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		if (Value < Ts[i])
			return i;
	return AlphaSize - 1;
	}

void DSS::Condense(const vector<float> &UnsortedValues, uint AlphaSize,
				   UNDEF_BINNING UB, uint BestDefaultLetter, uint &DefaultLetter,
				   float &MinValue, float &MedValue, float &MaxValue,
				   float &UndefFreq, vector<float> &BinTs)
	{
	BinTs.clear();
	vector<float> SortedValues;
	const uint N = SIZE(UnsortedValues);
	uint UndefCount = 0;
	for (uint i = 0; i < N; ++i)
		{
		float Value = UnsortedValues[i];
		if (Value == FLT_MAX)
			{
			++UndefCount;
			if (UB == UB_NeverUndefined)
				Die("Undefined");
			}
		else
			SortedValues.push_back(Value);
		}
	sort(SortedValues.begin(), SortedValues.end());
	const uint M = SIZE(SortedValues);
	MinValue = SortedValues[0];
	MedValue = SortedValues[M/2];
	MaxValue = SortedValues[M-1];
	UndefFreq = (float) UndefCount/M;

	const uint K = SIZE(SortedValues);
	asserta(K > 0);
	float FirstValue = SortedValues[0];
	float LastValue = SortedValues[K-1];
	asserta(FirstValue < LastValue);
	asserta(FirstValue != -FLT_MAX);
	asserta(LastValue != FLT_MAX);

	DefaultLetter = UINT_MAX;
	switch (UB)
		{
	case UB_IgnoreUndefined:
	case UB_NeverUndefined:
		DefaultLetter = UINT_MAX;
		break;

	case UB_UndefinedIsOnlyZero:
	case UB_UndefinedIsZeroOverload:
		DefaultLetter = 0;
		break;

	case UB_UndefinedIsDefaultLetter:
		asserta(BestDefaultLetter != UINT_MAX);
		DefaultLetter = BestDefaultLetter;
		break;

	default:
		asserta(false);
		}

	uint BinThresholdCount = GetBinThresholdCount(AlphaSize, UB);
	uint BinCount = BinThresholdCount + 1;
	for (uint i = 0; i < BinThresholdCount; ++i)
		{
		uint k = ((i+1)*K)/BinCount;
		float t = SortedValues[k];
		if (i > 0)
			{
			if (t == SortedValues[i-1])
				Warning("Condense tie for bin thresholds %u,%u at %.3g",
						i-1, i, t);
			asserta(t >= SortedValues[i-1]);
			}
		BinTs.push_back(t);
		}
	asserta(SIZE(BinTs) + 1 == BinCount);
	}
