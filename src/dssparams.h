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
	static bool m_InitDone;
	static vector<FEATURE> m_Features;
	static vector<float> m_Weights;

	static float ***m_ScoreMxs;
	static bool m_OwnScoreMxs;

	static float m_GapOpen;
	static float m_GapExt;
	static float m_FwdMatchScore;
	static float m_MinFwdScore;
	static float m_Omega;
	static float m_OmegaFwd;
	static string m_MKFPatternStr;
	static string m_MuPrefilterPatternStr;

	static int m_ParaMuGapOpen;
	static int m_ParaMuGapExt;

	static uint m_MKFL;
	static int m_MKF_X1;
	static int m_MKF_X2;
	static int m_MKF_MinHSPScore;
	static float m_MKF_MinMegaHSPScore;

public:
	static void Clear();
	static void SetDefaults();
	static void SetDefaults_Features();
	static void SetDefaults_Other();
	static void Init(DECIDE_MODE DM);

	static uint GetFeatureCount()
		{
		return SIZE(m_Features);
		}

	static void AddFeature(FEATURE F, double w)
		{
		m_Features.push_back(F);
		m_Weights.push_back(float(w));
		}

	static void LoadFeatures(const string &FN);
	static FEATURE LoadFeature(const string &FN);

	static void NormalizeWeights();
	static void AllocScoreMxs();
	static void SetScoreMxs();
	static void InitScoreMxs();
	static void ApplyWeights();
	};

uint GetPatternOnes(const string &Pattern);
