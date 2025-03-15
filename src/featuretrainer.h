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
	DSSParams m_Params;
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
	DSS m_DQ;
	DSS m_DR;

public:
	void SetFeature(FEATURE F);
	void SetAlphaSize(uint AS);
	void SetInput(const string &ChainsFN, const string &AlnsFN);
	void Train();
	void WriteSummary(FILE *f) const;
	void ToTsv(const string &FN) const;
	void FromTsv(const string &FN);

private:
	void SetLabelToChainIndex();
	void SetFloatValues();
	void SetJointFreqsPair(uint PairIndex);
	};
