#pragma once

#include "features.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "dss.h"
#include "logodds.h"

// FeatureTrainer and Trainer should be merged
class FeatureTrainer : public LogOdds
	{
public:
	FEATURE m_F = FEATURE(-1);
	const char *m_FeatureName = 0;
	vector<PDBChain *> m_Chains;
	bool m_IsInt = false;
	SeqDB m_Alns;
	int8_t m_MaxAbsi8 = 20;
	vector<float> m_FloatValues;
	uint m_UndefCount = 0;
	vector<uint> m_Letters;
	map<string, uint> m_LabelToChainIndex;
	vector<float> m_BinTs;
	DSS m_D;
	uint m_BestDefaultLetter = UINT_MAX;
	float m_BestDefaultValue = FLT_MAX;
	mutex m_Lock;
	uint m_Counter = UINT_MAX;
	vector<vector<float> > m_ScoreMx;

public:
	void SetFeature(FEATURE F, uint AlphaSize);
	void SetInput(const string &ChainsFN, const string &AlnsFN);
	void TrainLogOdds(bool IgnoreUndef);
	void WriteSummary(FILE *f) const;
	void ToTsv(const string &FN) const;
	void FromTsv(const string &FN);
	void ToSrc(FILE *f) const;
	void BinTsToSrc(FILE *f) const;
	void FreqsToSrc(FILE *f) const;
	void ScoreMxToSrc(FILE *f) const;
	void ScoreMxFromTsv(FILE *f);
	void TrainQuantization(bool UndefOverlap);
	float GetDefaultValue() const;

	void SetLabelToChainIndex();
	void SetFloatValues(bool IgnoreUndef);
	void UpdateJointCounts(uint PairIndex, bool IgnoreUndef);
	void AddValue(float Value);
	void SetUnalignedBackground(bool IgnoreUndef);
	void SetUnalignedBackgroundChain(const PDBChain &Chain,
		uint UndefLetter);
	float GetPctIdFromLabel(const string &Label) const;
	float GetAQFromLabel(const string &Label) const;
	};
