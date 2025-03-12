#include "myutils.h"
#include "dss.h"

/***
Undefined values:
	Letter is zero.
	DSS:GetFloat_F returns FLT_MAX.

Alphabet has AlphaSize (AS) letters.
Letter 0 is undefined.
That leaves AS-1 defined letters separated by AS-2 thresholds.
***/

uint DSS::ValueToInt(const vector<float> &Ts, float Value)
	{
	if (Value == FLT_MAX)
		return 0;
	const uint N = SIZE(Ts);
	for (uint i = 0; i < N; ++i)
		if (Value <= Ts[i])
			return i+1;
	return N+1;
	}

void DSS::Condense(const vector<float> &UnsortedValues, uint AlphaSize,
				   float &MinValue, float &MedValue, float &MaxValue, float &UndefFreq,
				   vector<float> &BinTs)
	{
	BinTs.clear();
	vector<float> SortedValues;
	const uint N = SIZE(UnsortedValues);
	uint UndefCount = 0;
	for (uint i = 0; i < N; ++i)
		{
		float Value = UnsortedValues[i];
		if (Value == FLT_MAX)
			++UndefCount;
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

	uint BinCount = AlphaSize - 1;
	for (uint i = 0; i + 1 < BinCount; ++i)
		{
		//uint k = ((i+1)*K)/(BinCount + 1);
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
	asserta(SIZE(BinTs) + 2 == AlphaSize);
	}
