#include "myutils.h"
#include "dss.h"

/***
Alphabet has AlphaSize (AS) letters.
Quantize by AS-1 thresholds
Undefined value can be FLT_MAX or in range.
***/

uint DSS::ValueToInt(float Value, uint AlphaSize, const vector<float> &Ts,
	uint DefaultLetter)
	{
	if (Value == FLT_MAX)
		{
		asserta(DefaultLetter < AlphaSize);
		return DefaultLetter;
		}

	asserta(SIZE(Ts) + 1 == AlphaSize);
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		if (Value < Ts[i])
			return i;
	return AlphaSize - 1;
	}

// UndefinedValue=FLT_MAX ignore undefineds
// Otherwise replace Value=FLT_MAX with UndefinedValue
void DSS::Quantize(const vector<float> &Values, uint AlphaSize,
	float UndefValue, vector<float> &BinTs)
	{
	BinTs.clear();

	const uint InputValueCount = SIZE(Values);
	vector<float> SortedValues;
	SortedValues.reserve(InputValueCount);
	uint UndefCount = 0;
	for (uint i = 0; i < InputValueCount; ++i)
		{
		float Value = Values[i];
		if (Value == FLT_MAX)
			Value = UndefValue;
		if (Value == FLT_MAX)
			continue;
		SortedValues.push_back(Value); 
		}
	sort(SortedValues.begin(), SortedValues.end());
	const uint K = SIZE(SortedValues);
	asserta(K > 0);
	float MinValue = SortedValues[0];
	float MaxValue = SortedValues[K-1];
	asserta(MinValue < MaxValue);
	asserta(MinValue != -FLT_MAX);
	asserta(MaxValue != FLT_MAX);

	for (uint i = 0; i + 2 < AlphaSize; ++i)
		{
		uint k = ((i+1)*K)/AlphaSize;
		float t = SortedValues[k];
		if (i > 0)
			{
			if (t == SortedValues[i-1])
				Die("Quantize tie for bin thresholds %u,%u at %.3g",
						i-1, i, t);
			asserta(t >= SortedValues[i-1]);
			}
		BinTs.push_back(t);
		}
	BinTs.push_back(MaxValue+1);
	asserta(SIZE(BinTs) + 1 == AlphaSize);
	}
