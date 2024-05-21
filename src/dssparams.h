#pragma once

#include "features.h"

class DSSParams
	{
private:
	vector<float> m_Weights;

public:
	string m_Desc;

	vector<FEATURE> m_Features;

	float m_GapOpen = FLT_MAX;
	float m_GapExt = FLT_MAX;
	float m_DALIw = FLT_MAX;
	float m_FwdMatchScore = FLT_MAX;
	float m_MinFwdScore = FLT_MAX;
	float m_Omega = FLT_MAX;
	uint m_MinU = 0;
	string m_PatternStr = "";
	float ***m_ScoreMxs = 0;
	bool m_USort = false;

	float m_EvalueSlope = -6.6f;
	float m_EvalueIntercept = 6.1f;
	uint m_Lambda = 32;

	float m_DBSize = FLT_MAX;
	uint m_MaxAccepts = UINT_MAX;
	uint m_MaxRejects = UINT_MAX;

// Not used
	uint m_MAXNQNR = INT_MAX;
	float m_X = FLT_MAX;
	float m_MAXFX = FLT_MAX;
	int m_MinDiagSpacing = INT_MAX;
	int m_MaxDiagSpacing = INT_MAX;
	float m_MinKmerScore = FLT_MAX;

public:
	void Clear()
		{
		m_Desc.clear();
		m_Features.clear();
		m_Weights.clear();

		m_GapOpen = FLT_MAX;
		m_GapExt = FLT_MAX;
		m_DALIw = FLT_MAX;
		m_FwdMatchScore = FLT_MAX;
		m_MinFwdScore = FLT_MAX;
		m_MAXNQNR = INT_MAX;
		m_MAXFX = FLT_MAX;
		m_MinDiagSpacing = INT_MAX;
		m_MaxDiagSpacing = INT_MAX;
		m_MinKmerScore = FLT_MAX;
		m_Omega = FLT_MAX;
		m_PatternStr = "";
		m_X = FLT_MAX;
		}

	void SetDefaults()
		{
		SetNamedParams("defaults");
		}

	void AddFeature(FEATURE F, double w)
		{
		m_Features.push_back(F);
		m_Weights.push_back(float(w));
		}

	void SetNamedParams(const string &Name);
	void NormalizeWeights();
	void WriteSummary(FILE *f) const;
	uint GetFeatureCount() const;
	void SetParam(const string &Name, float Value, bool AppendIfWeight);
	float GetParam(const string &Name) const;
	void SetFromCmdLine();
	uint GetFeatureIdx(FEATURE F) const;
	uint GetFeatureIdx_NoError(FEATURE F) const;
	void ToFev(FILE *f) const;
	void FromTsv(const string &FileName);
	void InitScoreMxs();
	void ApplyWeights();
	float ScoreToEvalue(float Score, uint QL) const;
	};

uint GetPatternOnes(const string &Pattern);
