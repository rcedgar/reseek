#pragma once

#include "features.h"

#define SLOPE_CALIB		0
#define GUMBEL_CALIB	0

const uint SCOP40_DBSIZE = 11211;

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
	string m_MuPrefPatternStr = "";
	float ***m_ScoreMxs = 0;
	bool m_OwnScoreMxs = false;

	float m_DBSize = 10000;

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
		m_MuPrefPatternStr = "";
		}

	void SetDefaults();

	void AddFeature(FEATURE F, double w)
		{
		m_Features.push_back(F);
		m_Weights.push_back(float(w));
		}

	void FromParamStr(const string &ParamStr);
	void NormalizeWeights();
	uint GetFeatureCount() const;
	void SetParam(const string &Name, float Value, bool AppendIfWeight);
	void SetIntParam(const string &Name, int Value);
	float GetParam(const string &Name) const;
	int GetIntParam(const string &Name) const;
	void SetDSSParams(DECIDE_MODE DM, uint DBSize);
	uint GetFeatureIdx(FEATURE F) const;
	uint GetFeatureIdx_NoError(FEATURE F) const;
	void ToFev(FILE *f, bool nl) const;
	void FromTsv(const string &FileName);
	void InitScoreMxs();
	void ApplyWeights();
	float GetEvalue(float TS) const;
	};

uint GetPatternOnes(const string &Pattern);
