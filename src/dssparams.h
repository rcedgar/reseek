#pragma once

#include "features.h"

#define SLOPE_CALIB		0
#define GUMBEL_CALIB	0

enum ALGO_MODE
	{
	AM_Invalid,
	AM_Fast,
	AM_Sensitive,
	AM_VerySensitive
	};

enum DECIDE_MODE
	{
	DM_Invalid,
	DM_AlwaysFast,
	DM_AlwaysSensitive,
	DM_AlwaysVerysensitive,
	DM_DefaultFast,
	DM_DefaultSensitive,
	DM_UseCommandLineOption
	};

class DSSParams
	{
public:
	vector<FEATURE> m_Features;
	vector<float> m_Weights;

	float m_GapOpen = FLT_MAX;
	float m_GapExt = FLT_MAX;
	float m_FwdMatchScore = FLT_MAX;
	float m_MinFwdScore = FLT_MAX;
	float m_Omega = FLT_MAX;
	float m_OmegaFwd = FLT_MAX;
	string m_MKFPatternStr = "";
	string m_MuPrefilterPatternStr = "";
	float ***m_ScoreMxs = 0;
	bool m_OwnScoreMxs = false;

	bool m_UsePara = true;
	int m_ParaMuGapOpen = 2;
	int m_ParaMuGapExt = 1;

	float m_Evalue_old_linear_Slope = -6.6f;
	float m_Evalue_linear_Intercept = 6.1f;

	float m_Evalue_Gumbel_mu = 2.5f;
	float m_Evalue_Gumbel_beta = 0.613f;

	float m_Evalue_linear_m = 20.5f;
	float m_Evalue_linear_b = 2.9f;

	float m_Evalue_a = 4.0f;
	float m_Evalue_b = -43.0f;

	bool m_AAOnly = false;

	uint m_MKFL = UINT_MAX;
	int m_MKF_X1 = INT_MAX;
	int m_MKF_X2 = INT_MAX;
	int m_MKF_MinHSPScore = INT_MAX;
	float m_MKF_MinMegaHSPScore = FLT_MAX;

public:
	static const uint m_MuFeatureCount = 3;
	static const FEATURE m_MuFeatures[m_MuFeatureCount];
	static const uint m_MuAlphaSizes[3];
	static const uint m_MuAlphaSize;

public:
	~DSSParams();

public:
	void Clear()
		{
		m_Features.clear();
		m_Weights.clear();

		m_GapOpen = FLT_MAX;
		m_GapExt = FLT_MAX;
		m_FwdMatchScore = FLT_MAX;
		m_MinFwdScore = FLT_MAX;
		m_Omega = FLT_MAX;
		m_MKFPatternStr = "";
		m_MuPrefilterPatternStr = "";
		}

	void SetDefaults();
	void SetDefaults_Features();
	void SetDefaults_Other();

	void AddFeature(FEATURE F, double w)
		{
		m_Features.push_back(F);
		m_Weights.push_back(float(w));
		}

	void LoadFeatures();
	FEATURE LoadFeature(const string &FN);

	void FromParamStr(const string &ParamStr);
	void NormalizeWeights();
	uint GetFeatureCount() const;
	void SetParam(const string &Name, float Value, bool AppendIfWeight);
	void SetIntParam(const string &Name, int Value);
	float GetParam(const string &Name) const;
	int GetIntParam(const string &Name) const;
	void SetDSSParams(DECIDE_MODE DM);
	uint GetFeatureIdx(FEATURE F) const;
	uint GetFeatureIdx_NoError(FEATURE F) const;
	void ToFev(FILE *f, bool nl) const;
	void FromTsv(const string &FileName);
	void AllocScoreMxs();
	void SetScoreMxs();
	void InitScoreMxs();
	void ApplyWeights();
	//float GetEvalue(float TS) const;
	//float GetEvalueSlope(float TS, float m, float b) const;
	//float GetEvalueGumbel(float TS, float mu, float beta) const;
	//float GetEvalueOldLinear(float TS) const;
	};

uint GetPatternOnes(const string &Pattern);
