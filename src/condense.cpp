#include "myutils.h"
#include "dss.h"

/***
Alphabet has AlphaSize (AS) letters.
Quantize by AS-1 thresholds
Undefined value is FLT_MAX

Dedicated undefined letter
--------------------------
Defined letters are 0..AS-2, undefined is AS-1

AS-2 thresholds in BinTs
BinTs[AS-2] = MaxValue+1
	MaxValue => AS-2

Overloaded letter for undefined
-------------------------------
Defined letters are 0..AS-1, any one letter is overloaded

AS-1 thresholds in BinTs
BinTs[AS-1] = MaxValue+1
	MaxValue => AS-1
***/

uint DSS::ValueToInt(float Value, uint AlphaSize, const vector<float> &Ts,
	uint DefaultLetter)
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

void DSS::Quantize(const vector<float> &Values, uint AlphaSize,
	bool OverloadUndefined, uint UndefinedLetter, vector<float> &BinTs)
	{
	asserta(UndefinedLetter < AlphaSize);
	BinTs.clear();

	const uint N = SIZE(Values);
	vector<float> SortedValues;
	SortedValues.reserve(N);
	uint UndefCount = 0;
	for (uint i = 0; i < N; ++i)
		{
		float Value = Values[i];
		if (Value != FLT_MAX)
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
