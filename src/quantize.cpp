#include "myutils.h"
#include "featuretrainer.h"

void FeatureTrainer::Quantize(const vector<float> &SortedValues,
	uint AlphaSize, vector<float> &BinTs)
	{
	BinTs.clear();

	const uint K = SIZE(SortedValues);
	asserta(K > 0);
	float MinValue = SortedValues[0];
	float MaxValue = SortedValues[K-1];
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		{
		uint k = ((i+1)*K)/AlphaSize;
		float t = SortedValues[k];
		if (i > 0)
			{
			if (feq(t, BinTs[i-1]))
				Log("Quantize tie for bin thresholds %u,%u at %.3g\n",
						i-1, i, t);
			QuantizeUniques(SortedValues, AlphaSize, BinTs);
			asserta(SIZE(BinTs) + 1 == AlphaSize);
			return;
			}
		BinTs.push_back(t);
		}
	asserta(SIZE(BinTs) + 1 == AlphaSize);
	}

void FeatureTrainer::QuantizeUniques(const vector<float> &SortedValues,
	uint AlphaSize, vector<float> &BinTs)
	{
	BinTs.clear();
	vector<float> UniqueFloatValues;
	vector<uint> UniqueFloatCounts;
	
	const uint N = SIZE(SortedValues);
	asserta(N > 100);
	float UniqueValue = SortedValues[0];
	uint Count = 1;
	for (uint i = 1; i < N; ++i)
		{
		float Value = SortedValues[i];
		if (Value == UniqueValue)
			++Count;
		else
			{
			asserta(Value > UniqueValue);
			UniqueFloatValues.push_back(UniqueValue);
			UniqueFloatCounts.push_back(Count);
			UniqueValue = Value;
			Count = 1;
			}
		}
	UniqueFloatValues.push_back(UniqueValue);
	UniqueFloatCounts.push_back(Count);
	uint UniqueValueCount = SIZE(UniqueFloatValues);
	asserta(SIZE(UniqueFloatCounts) == UniqueValueCount);
	asserta(UniqueValueCount >= AlphaSize);
	uint TargetCountPerBin = uint(double(N)/AlphaSize + 0.5);
	vector<float> TmpValues;
	for (uint i = 0; i < UniqueValueCount; ++i)
		{
		float Value = UniqueFloatValues[i];
		uint n = min(UniqueFloatCounts[i], TargetCountPerBin/2);
		for (uint j = 0; j < n; ++j)
			TmpValues.push_back(Value);
		}

	const uint K = SIZE(TmpValues);
	asserta(K > 0);
	float MinValue = TmpValues[0];
	float MaxValue = TmpValues[K-1];
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		{
		uint k = ((i+1)*K)/AlphaSize;
		float t = TmpValues[k];
		if (i > 0)
			{
			if (feq(t, BinTs[i-1]))
				Die("QuantizeUniques tie for bin thresholds %u,%u at %.3g\n",
						i-1, i, t);
			}
		BinTs.push_back(t);
		}
	asserta(SIZE(BinTs) + 1 == AlphaSize);
	}
