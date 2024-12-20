#pragma once

#include "features.h"

#define SLOPE_CALIB		0
#define GUMBEL_CALIB	0

class DSSParams
	{
public:
	string m_Desc;

	vector<FEATURE> m_Features;
	vector<float> m_Weights;

	float m_GapOpen = FLT_MAX;
	float m_GapExt = FLT_MAX;
	float m_DALIw = FLT_MAX;
	float m_FwdMatchScore = FLT_MAX;
	float m_MinFwdScore = FLT_MAX;
	float m_MinComboFwdScore = FLT_MAX;
	float m_Omega = FLT_MAX;
	float m_OmegaFwd = FLT_MAX;
	uint m_MinU = 0;
	string m_PatternStr = "";
	float ***m_ScoreMxs = 0;
	bool m_USort = false;
	//bool m_ComboScoreOnly = false;
	//bool m_UseComboPath = false;

	uint m_Lambda = 32;

	float m_DBSize = 10000;
	uint m_MaxAccepts = UINT_MAX;
	uint m_MaxRejects = UINT_MAX;

	bool m_UsePara = true;
	//int m_ParaComboGapOpen = 5;
	//int m_ParaComboGapExt = 1;
	int m_ParaComboGapOpen = 2;
	int m_ParaComboGapExt = 1;

	float m_Evalue_old_linear_Slope = -6.6f;
	float m_Evalue_linear_Intercept = 6.1f;

	float m_Evalue_Gumbel_mu = 2.5f;
	float m_Evalue_Gumbel_beta = 0.613f;

// Superfamily: Linear fit to -log(P) m=20.5 b=2.89
// Fold:        Linear fit to -log(P) m=26.6 b=2.42
	float m_Evalue_linear_m = 20.5f;
	float m_Evalue_linear_b = 2.9f;

	float m_Evalue_a = 4.0f;
	float m_Evalue_b = -43.0f;

	bool m_AAOnly = false;

public:
	static vector<FEATURE> m_ComboFeatures;
	static vector<uint> m_ComboAlphaSizes;
	static uint m_ComboAlphaSize;

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
		m_Omega = FLT_MAX;
		m_PatternStr = "";
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
	void FromParamStr(const string &ParamStr);
	void NormalizeWeights();
	void WriteSummary(FILE *f) const;
	uint GetFeatureCount() const;
	void SetParam(const string &Name, float Value, bool AppendIfWeight);
	void SetIntParam(const string &Name, int Value);
	float GetParam(const string &Name) const;
	int GetIntParam(const string &Name) const;
	void SetFromCmdLine(uint DBSize);
	uint GetFeatureIdx(FEATURE F) const;
	uint GetFeatureIdx_NoError(FEATURE F) const;
	void ToFev(FILE *f, bool nl) const;
	void FromTsv(const string &FileName);
	void InitScoreMxs();
	void ApplyWeights();
	float GetEvalue(float TS) const;
	float GetEvalueSlope(float TS, float m, float b) const;
	float GetEvalueGumbel(float TS, float mu, float beta) const;
	float GetEvalueOldLinear(float TS) const;

public:
	static void SetComboFeatures(const vector<FEATURE> &Fs);
	};

uint GetPatternOnes(const string &Pattern);
