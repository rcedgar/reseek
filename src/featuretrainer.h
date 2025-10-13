#pragma once

#include "features.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "dss.h"
#include "logodds.h"
#include "peaker.h"

// FeatureTrainer and Trainer should be merged
class FeatureTrainer : public LogOdds
	{
public:
	FEATURE m_F = FEATURE(-1);
	const char *m_FeatureName = 0;
	bool m_IsInt = false;
	string m_UndefStyle = "ERROR";
	vector<PDBChain *> m_Chains;
	vector<vector<float> > m_ChainFloatSeqVec;
	vector<vector<uint> > m_ChainLetterSeqVec;
	map<string, uint> m_LabelToChainIndex;

	vector<string> m_AlnLabels;
	vector<string> m_AlnRows;
	vector<uint> m_AlnChainIdxs;
	vector<bool> m_AlnTPs;
	vector<uint> m_AlnOpenCounts;
	vector<uint> m_AlnExtCounts;
	uint m_TPCount = 0;
	uint m_FPCount = 0;

	vector<float> m_AlnSubstScores;
	vector<float> m_AlnScores;

	vector<float> m_ROCStepScores;
	vector<float> m_ROCStepTPfs;
	vector<float> m_ROCStepFPfs;
	vector<uint> m_ROCOrder;

	DSS m_D;
	vector<float> m_SortedFloatValues;
	vector<uint> m_Letters;
	uint m_UndefCount = 0;

	int8_t m_MaxAbsi8 = 20;
	vector<float> m_BinTs;
	uint m_BestDefaultLetter = UINT_MAX;
	float m_BestDefaultValue = FLT_MAX;
	vector<vector<float> > m_ScoreMx;

	float m_OpenPenalty = FLT_MAX;
	float m_ExtPenalty = FLT_MAX;
	float m_Area = FLT_MAX;

	Peaker m_Peaker;

	mutex m_Lock;

public:
	void ReadChains(const string &ChainsFN);
	void ReadAlns(const string &AlnsFN, bool TPs);
	void AddAln(const string &Label1, const string &Row1,
		const string &Label2, const string &Row2, bool TP);
	void SetFeature(FEATURE F, uint AlphaSize);
	void TrainLogOdds(bool IgnoreUndef);
	void WriteSummary(FILE *f) const;
	void ToTsv(const string &FN) const;
	void FromTsv(const string &FN);
	void ToSrc(FILE *f) const;
	void BinTsToSrc(FILE *f) const;
	void FreqsToSrc(FILE *f) const;
	void ScoreMxToSrc(FILE *f) const;
	void ScoreMxFromTsv(FILE *f);
	uint GetUngappedLength(const string &Row) const;
	void TruncLabel(string &Label);
	
	void TrainFloat_UndefOverlap();
	void TrainFloat_UndefDistinct();
	void TrainInt_UndefOverlap();
	void TrainInt_UndefDistinct();
	float GetDefaultValue() const;
	float GetMaxDefinedValue() const;
		
	void SetLabelToChainIndex();
	void SetFloatValues(bool IgnoreUndef, float ReplaceUndefValue);
	void UpdateJointCounts(uint PairIndex, bool IgnoreUndef);
	void SetUnalignedBackground(bool IgnoreUndef);
	void SetUnalignedBackgroundChain(const PDBChain &Chain,
		uint UndefLetter);
	void SetChainFloatSeqs(float ReplaceUndefValue);
	void SetChainLetterSeqs();
	void SetChainLetterSeqs_Float();
	void SetChainLetterSeqs_Int();
	void GetGapCounts(const string &Row1, const string &Row2,
		uint &Opens, uint &Exts) const;
	uint GetBestUndefLetter() const;

	void SetAlnSubstScores();
	float GetAlnSubstScore(uint AlnIdx);
	void SetAlnScoresAndArea(float OpenPenalty, float ExtPenalty);

	uint GetAlnCount() const { return SIZE(m_AlnTPs); }

	void LogROCStepsAndArea();
	void SetArea();

	void OptimizeGapPenalties();

public:
	void Quantize(const vector<float> &Values, uint AlphaSize,
		vector<float> &BinTs);
	void QuantizeUniques(const vector<float> &SortedValues,
		uint AlphaSize, vector<float> &BinTs);
	
	static double EvalArea(const vector<double> &xv);
	};
