#pragma once

#include "features.h"

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
private:
	DSSParams() { Die("DSSParams()"); }
	~DSSParams() { Die("~DSSParams()"); }

public:
	static vector<FEATURE> m_Features;
	static vector<float> m_Weights;
	static float ***m_ScoreMxs;

	static float m_GapOpen;
	static float m_GapExt;
	static float m_MinFwdScore;
	static float m_Omega;
	static float m_OmegaFwd;
	static string m_MKFPatternStr;

	static int m_ParaMuGapOpen;
	static int m_ParaMuGapExt;

	static uint m_rsb_size;
	static int m_PrefilterMinKmerPairScore;
	static uint m_MKFL;
	static int m_MKF_X1;
	static int m_MKF_X2;
	static int m_MKF_MinMuHSPScore;
	static float m_MKF_MinMegaHSPScore;
	static float m_dpw;
	static float m_lddtw;
	static float m_ladd;
	static float m_revtsw;

	static int m_Density_W;
	static int m_Density_w;
	static int m_SSDensity_W;
	static int m_SSDensity_w;
	static float m_Density_Radius;
	static float m_NU_ND_Radius;
	static int m_NEN_W;
	static int m_NEN_w;
	static int m_NUDX_W;
	static float m_DefaultNENDist;
	static float m_SSDensity_epsilon;
	static uint m_SSE_MinLength;
	static uint m_SSE_Margin;
	static uint m_PMDelta;

public:
	static void Init(DECIDE_MODE DM);
	static void SetStandardFeatures();
	static void SetDefaultNonFeatureTunableParams();
	//static void SetTunableParamFromStr(const string &Name, const string &Value);
	static void SetParamsFromStr(const string &Str);
	static void SetAlgoMode(DECIDE_MODE DM);
	static void GetParamsStr(string &Str);
	static void LogMe();

	static uint GetFeatureCount()
		{
		return SIZE(m_Features);
		}

	static void AddFeature(FEATURE F, float w)
		{
		m_Features.push_back(F);
		m_Weights.push_back(w);
		}

	static void DumpFeature(FEATURE F, const string &FN);

	static void LoadFeatures(const string &FN);
	static FEATURE LoadFeature(const string &FN, float Weight);

	static void NormalizeWeights();
	static void AllocScoreMxs();
	static void FreeScoreMxs();
	static void SetScoreMxsFromFeatures();
	static void CreateWeightedScoreMxs();
	static void UpdateWeights(const vector<float> &Weights);

	static uint GetAlphaSize(FEATURE F, bool FailOk = false);
	static void GetBinTs(FEATURE F, vector<float> &Bins);
	static uint ValueToInt_Feature(FEATURE F, float Value);

	static void OverwriteUnweightedScoreMx(FEATURE F,
		vector<vector<float> > &ScoreMx);
	static void OverwriteBinTs(FEATURE F,
		const vector<float> &BinTs);

	static void InitParamData();
	static void LogParamData();
	static bool ParamIsDefault(const string &Name);
	static const char *GetParamDefaultStr(const string &Name, string &DefaultStr);
	static const char *GetParamValueStr(const string &Name, string &ValueStr);
	static int GetIntParam(const string &Name);
	static uint GetUintParam(const string &Name);
	static float GetFloatParam(const string &Name);
	static void SetParam(const string &Name, const string &Value);
	};

uint GetPatternOnes(const string &Pattern);
