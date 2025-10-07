#pragma once

#include "features.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "dss.h"
#include "logodds.h"
#include "undef_binning.h"

// FeatureTrainer and Trainer should be merged
class FeatureTrainer : public LogOdds
	{
public:
	FEATURE m_F = FEATURE(-1);
	const char *m_FeatureName = 0;
	float m_MinAQ = 0;
	float m_MaxAQ = 1;
	float m_MinPctId = 0;
	float m_MaxPctId = 100;
	vector<PDBChain *> m_Chains;
	bool m_IsInt = false;
	SeqDB m_Alns;
	int8_t m_MaxAbsi8 = 20;
	vector<float> m_FloatValues;
	vector<uint> m_Letters;
	float m_MinValue = FLT_MAX;
	float m_MedValue = -FLT_MAX;
	float m_MaxValue = -FLT_MAX;
	float m_UndefFreq = -FLT_MAX;
	map<string, uint> m_LabelToChainIndex;
	vector<float> m_BinTs;
	DSS m_D;
	//DSS m_DQ;
	//DSS m_DR;
	UNDEF_BINNING m_UB = UB_Invalid;
	uint m_BestDefaultLetter = UINT_MAX;
	uint m_ExcludedPairCount = UINT_MAX;
	mutex m_Lock;
	uint m_Counter = UINT_MAX;
	vector<vector<float> > m_ScoreMx;

public:
	void SetFeature(FEATURE F);
	void SetOptionsFromCmdLine();
	void SetAlphaSize(uint AS, UNDEF_BINNING UB, uint DefaultLetter);
	void SetInput(const string &ChainsFN, const string &AlnsFN);
	void Train();
	void WriteSummary(FILE *f) const;
	void ToTsv(const string &FN) const;
	void FromTsv(const string &FN);
	void ToSrc(FILE *f) const;
	void BinTsToSrc(FILE *f) const;
	void FreqsToSrc(FILE *f) const;
	void ScoreMxToSrc(FILE *f) const;
	void ScoreMxFromTsv(FILE *f);

private:
	void SetLabelToChainIndex();
	void SetFloatValues();
	void UpdateJointCounts(uint PairIndex);
	bool IncludePair(const string &QLabel) const;
	void AddValue(float Value);
	void SetUnalignedBackground();
	void SetUnalignedBackgroundChain(const PDBChain &Chain);
	float GetPctIdFromLabel(const string &Label) const;
	float GetAQFromLabel(const string &Label) const;
	};
