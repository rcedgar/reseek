#include "myutils.h"
#include "featuretrainer2.h"

void FeatureTrainer2::Quantize(
	const vector<float> &InputSortedValues,
	vector<float> &BinTs,
	float &UndefReplaceValue)
	{
	BinTs.clear();
	const uint InputCount = SIZE(InputSortedValues);
	vector<float> SortedValues;
	SortedValues.reserve(InputCount);

	switch (m_QS)
		{
	case QS_DiscardUndef:
		Quantize_DiscardUndef(InputSortedValues, BinTs);
		break;

	case QS_UndefDistinct:
		Quantize_UndefDistinct(InputSortedValues, BinTs);
		break;

	case QS_UndefNotSpecialCase:
		Quantize_UndefNotSpecialCase(InputSortedValues, BinTs);
		break;

	case QS_UndefReplaceUser:
		Quantize_UndefReplaceUser(InputSortedValues,
			BinTs, UndefReplaceValue);
		break;

	default:
		asserta(false);
		}

	bool Tie = false;
	for (uint i = 1; i + 1 < m_AlphaSize; ++i)
		{
		if (feq(BinTs[i-1], BinTs[i]))
			Warning("Quantize QS=%d tie bin %u T=%.3g,%.3g\n",
				int(m_QS), i, BinTs[i-1], BinTs[i]);
		}
	}

void FeatureTrainer2::Quantize_DiscardUndef(
	const vector<float> &InputSortedValues,
	vector<float> &BinTs)
	{
	BinTs.clear();
	const uint InputCount = SIZE(InputSortedValues);
	vector<float> SortedValues;
	SortedValues.reserve(InputCount);

	for (uint i = 0; i < InputCount; ++i)
		{
		float Value = InputSortedValues[i];
		if (Value != FLT_MAX)
			SortedValues.push_back(Value);
		}

	const uint K = SIZE(SortedValues);
	asserta(K > 0);
	for (uint i = 0; i + 1 < m_AlphaSize; ++i)
		{
		uint k = ((i+1)*K)/m_AlphaSize;
		float t = SortedValues[k];
		BinTs.push_back(t);
		}
	asserta(SIZE(BinTs) + 1 == m_AlphaSize);
	}

void FeatureTrainer2::Quantize_UndefNotSpecialCase(
	const vector<float> &InputSortedValues,
	vector<float> &BinTs)
	{
	BinTs.clear();
	const uint InputCount = SIZE(InputSortedValues);
	const vector<float> &SortedValues = InputSortedValues;
	const uint K = SIZE(SortedValues);
	asserta(K > 0);
	asserta(K == InputCount);
	for (uint i = 0; i + 1 < m_AlphaSize; ++i)
		{
		uint k = ((i+1)*K)/m_AlphaSize;
		float t = SortedValues[k];
		BinTs.push_back(t);
		}
	asserta(SIZE(BinTs) + 1 == m_AlphaSize);
	}

void FeatureTrainer2::Quantize_UndefDistinct(
	const vector<float> &InputSortedValues,
	vector<float> &BinTs)
	{
	BinTs.clear();
	const uint InputCount = SIZE(InputSortedValues);
	vector<float> SortedValues;
	SortedValues.reserve(InputCount);

// Discard undefs, they are special-cased
	for (uint i = 0; i < InputCount; ++i)
		{
		float Value = InputSortedValues[i];
		if (Value != FLT_MAX)
			SortedValues.push_back(Value);
		}

	const uint K = SIZE(SortedValues);
	asserta(K > 0);
// Create m_AlphaSize-1 bins (AS-2 thresholds) for defined letters
	for (uint i = 0; i + 2 < m_AlphaSize; ++i)
		{
		uint k = ((i+1)*K)/(m_AlphaSize - 1);
		float t = SortedValues[k];
		BinTs.push_back(t);
		}
	asserta(SIZE(BinTs) + 2 == m_AlphaSize);

// Special case bin threshold for undef, Letter=m_AlphaSize-1
// All values except FLT_MAX will be < this threshold
// *MUST* use < *NOT* <= in ValueToInt_xxx
	BinTs.push_back(FLT_MAX);
	}

void FeatureTrainer2::Quantize_UndefReplaceUser(
	const vector<float> &InputSortedValues,
	vector<float> &BinTs,
	float ReplaceValue)
	{
	BinTs.clear();
	const uint InputCount = SIZE(InputSortedValues);

	vector<float> SortedValues;
	SortedValues.reserve(InputCount);
	for (uint i = 0; i < InputCount; ++i)
		{
		float Value = InputSortedValues[i];
		if (Value == FLT_MAX)
			SortedValues.push_back(ReplaceValue);
		else
			SortedValues.push_back(Value);
		}
	sort(SortedValues.begin(), SortedValues.end());

	const uint K = SIZE(SortedValues);
	asserta(K == InputCount);
	asserta(K > 0);
	for (uint i = 0; i + 1 < m_AlphaSize; ++i)
		{
		uint k = ((i+1)*K)/m_AlphaSize;
		float t = SortedValues[k];
		BinTs.push_back(t);
		}
	asserta(SIZE(BinTs) + 1 == m_AlphaSize);
	}
