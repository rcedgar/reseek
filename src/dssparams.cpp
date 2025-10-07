#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "sort.h"
#include "prefiltermuparams.h"

bool DSSParams::m_InitDone = false;
vector<FEATURE> DSSParams::m_Features;
vector<float> DSSParams::m_Weights;

float DSSParams::m_GapOpen = FLT_MAX;
float DSSParams::m_GapExt = FLT_MAX;
float DSSParams::m_FwdMatchScore = FLT_MAX;
float DSSParams::m_MinFwdScore = FLT_MAX;
float DSSParams::m_Omega = FLT_MAX;
float DSSParams::m_OmegaFwd = FLT_MAX;
string DSSParams::m_MKFPatternStr = "";
string DSSParams::m_MuPrefilterPatternStr = "";
float ***DSSParams::m_ScoreMxs = 0;
bool DSSParams::m_OwnScoreMxs = false;

int DSSParams::m_ParaMuGapOpen = 2;
int DSSParams::m_ParaMuGapExt = 1;

uint DSSParams::m_MKFL = UINT_MAX;
int DSSParams::m_MKF_X1 = INT_MAX;
int DSSParams::m_MKF_X2 = INT_MAX;
int DSSParams::m_MKF_MinHSPScore = INT_MAX;
float DSSParams::m_MKF_MinMegaHSPScore = FLT_MAX;

static ALGO_MODE GetAlgoModeFromCommandLine(ALGO_MODE DefaultMode)
	{
	if (optset_fast)
		return AM_Fast;
	else if (optset_sensitive)
		return AM_Sensitive;
	else if (optset_verysensitive)
		return AM_VerySensitive;
	if (DefaultMode == AM_Invalid)
		Die("Must set -fast, -sensitive or -verysensitive");
	return DefaultMode;
	}

static ALGO_MODE GetAlgoMode(DECIDE_MODE DM)
	{
	switch (DM)
		{
	case DM_AlwaysFast:				return AM_Fast;
	case DM_AlwaysSensitive:		return AM_Sensitive;
	case DM_AlwaysVerysensitive:	return AM_VerySensitive;
	case DM_UseCommandLineOption:	return GetAlgoModeFromCommandLine(AM_Invalid);
	case DM_DefaultFast:			return GetAlgoModeFromCommandLine(AM_Fast);
	case DM_DefaultSensitive:		return GetAlgoModeFromCommandLine(AM_Sensitive);
		}
	asserta(false);
	return AM_Invalid;
	}

void DSSParams::Init(DECIDE_MODE DM)
	{
	asserta(!m_InitDone);
	SetDefaults();

	ALGO_MODE AM = GetAlgoMode(DM);

	switch (AM)
		{
	case AM_Fast:
		m_Omega = 22;
		m_OmegaFwd = 50;
		m_MKFL = 500;
		m_MKF_X1 = 8;
		m_MKF_X2 = 8;
		m_MKF_MinHSPScore = 50;
		m_MKF_MinMegaHSPScore = -4;
		break;

	case AM_Sensitive:
		m_Omega = 12;
		m_OmegaFwd = 20;
		m_MKFL = 600;
		m_MKF_X1 = 8;
		m_MKF_X2 = 8;
		m_MKF_MinHSPScore = 50;
		m_MKF_MinMegaHSPScore = -4;
		break;

	case AM_VerySensitive:
		m_Omega = 0;
		m_OmegaFwd = 0;
		m_MKFL = 99999;
		m_MKF_X1 = 99999;
		m_MKF_X2 = 99999;
		m_MKF_MinHSPScore = 0;
		m_MKF_MinMegaHSPScore = -99999;
		m_MinFwdScore = 0;
		break;

	default:
		asserta(false);
		}

	m_MKFPatternStr = "111";
	m_MuPrefilterPatternStr = string(prefiltermu_pattern);

	const int MINUS = -1; // for visual emphasis here
	if (optset_omega) { m_Omega = (float) opt(omega);  }
	if (optset_omegafwd) { m_OmegaFwd = (float) opt(omegafwd); }
	if (optset_minfwdscore) { m_MinFwdScore = float(opt(minfwdscore)); }
	if (optset_gapopen) { m_GapOpen =  MINUS*float(opt(gapopen)); }
	if (optset_gapopen) { m_GapExt = MINUS*float(opt(gapext)); }
	if (optset_para_mugapopen) { m_ParaMuGapOpen = opt(para_mugapopen); }
	if (optset_para_mugapext) { m_ParaMuGapExt = opt(para_mugapext); }
	if (optset_minhsp) { m_MKF_MinHSPScore = opt(minhsp); }
	if (optset_minmegahsp) { m_MKF_MinMegaHSPScore = float(opt(minmegahsp)); }
	if (optset_xdrop1) { m_MKF_X1 = int(opt(xdrop1)); }
	if (optset_xdrop2) { m_MKF_X2 = int(opt(xdrop2)); }
	if (optset_mkfl) { m_MKFL = int(opt(mkfl)); }

	if (m_GapOpen > 0 || m_GapExt > 0)
		Die("open=%.3g ext=%.3g, gap penalties must be >= 0",
		  opt(gapopen), opt(gapext));

	InitScoreMxs();
	m_InitDone = true;
	}

void DSSParams::NormalizeWeights()
	{
	float Sum = 0;
	const uint N = SIZE(m_Weights);
	for (uint Idx = 0; Idx < N; ++Idx)
		Sum += m_Weights[Idx];
	float Sum2 = 0;
	for (uint Idx = 0; Idx < N; ++Idx)
		{
		float w = m_Weights[Idx]/Sum;
		m_Weights[Idx] = w;
		Sum2 += w;
		}
	asserta(feq(Sum2, 1.0f));
	}

void DSSParams::InitScoreMxs()
	{
	if (m_ScoreMxs != 0)
		return;
	uint FeatureCount = GetFeatureCount();
	m_ScoreMxs = myalloc(float **, FEATURE_COUNT);
	for (uint i = 0; i < FEATURE_COUNT; ++i)
		m_ScoreMxs[i] = 0;
	for (uint Idx = 0; Idx < FeatureCount; ++Idx)
		{
		FEATURE F = m_Features[Idx];
		asserta(uint(F) < FEATURE_COUNT);
		uint AS = g_AlphaSizes2[F];
		asserta(m_ScoreMxs[F] == 0);
		m_ScoreMxs[F] = myalloc(float *, AS);
		for (uint Letter1 = 0; Letter1 < AS; ++Letter1)
			{
			m_ScoreMxs[F][Letter1] = myalloc(float, AS);
#if DEBUG
			for (uint Letter2 = 0; Letter2 < AS; ++Letter2)
				m_ScoreMxs[F][Letter1][Letter2] = FLT_MAX;
#endif
			}
		}
	ApplyWeights();
	m_OwnScoreMxs = true;
	}

void DSSParams::SetScoreMxs()
	{
	AllocScoreMxs();
	ApplyWeights();
	}

void DSSParams::AllocScoreMxs()
	{
	asserta(m_ScoreMxs == 0);
	uint FeatureCount = GetFeatureCount();
	m_ScoreMxs = myalloc(float **, FEATURE_COUNT);
	for (uint i = 0; i < FEATURE_COUNT; ++i)
		m_ScoreMxs[i] = 0;

	for (uint Idx = 0; Idx < FeatureCount; ++Idx)
		{
		FEATURE F = m_Features[Idx];
		asserta(uint(F) < FEATURE_COUNT);
		uint AS = DSS::GetAlphaSize(F);
		asserta(m_ScoreMxs[F] == 0);
		m_ScoreMxs[F] = myalloc(float *, AS);
		for (uint Letter1 = 0; Letter1 < AS; ++Letter1)
			{
			m_ScoreMxs[F][Letter1] = myalloc(float, AS);
#if DEBUG
			for (uint Letter2 = 0; Letter2 < AS; ++Letter2)
				m_ScoreMxs[F][Letter1][Letter2] = FLT_MAX;
#endif
			}
		}
	ApplyWeights();
	m_OwnScoreMxs = true;
	}

void DSSParams::ApplyWeights()
	{
	asserta(m_ScoreMxs != 0);
	uint FeatureCount = GetFeatureCount();
	for (uint Idx = 0; Idx < FeatureCount; ++Idx)
		{
		FEATURE F = m_Features[Idx];
		asserta(uint(F) < FEATURE_COUNT);
		float w = m_Weights[Idx];
		uint AS = g_AlphaSizes2[F];
		if (AS == 0)
			Die("Feature %s not supported", FeatureToStr(F));
		m_ScoreMxs[F] = myalloc(float *, AS);
		for (uint Letter1 = 0; Letter1 < AS; ++Letter1)
			{
			m_ScoreMxs[F][Letter1] = myalloc(float, AS);
			for (uint Letter2 = 0; Letter2 < AS; ++Letter2)
				m_ScoreMxs[F][Letter1][Letter2] = w*g_ScoreMxs2[F][Letter1][Letter2];
			}
		}
	}

void DSSParams::SetDefaults()
	{
	AddFeature(FEATURE_AA,			0.398145);
	AddFeature(FEATURE_NENDist,		0.129367);
	AddFeature(FEATURE_Conf,		0.202354);
	AddFeature(FEATURE_NENConf,		0.149383);
	AddFeature(FEATURE_RENDist,	0.0937677);
	AddFeature(FEATURE_DstNxtHlx,	0.00475462);
	AddFeature(FEATURE_StrandDens,	0.0183853);
	AddFeature(FEATURE_NormDens,	0.00384384);

	m_GapOpen = -0.685533f;
	m_GapExt = -0.051881f;
	m_FwdMatchScore = 0.1f;
	m_MinFwdScore = 7.0f;
	m_Omega = 29;
	m_OmegaFwd = 29;
	m_MuPrefilterPatternStr = "1110011";
	m_MKFPatternStr = "111";
	}

void DSSParams::Clear()
	{
	m_Features.clear();
	m_Weights.clear();
	m_GapOpen = FLT_MAX;
	m_GapExt = FLT_MAX;
	m_FwdMatchScore = FLT_MAX;
	m_MinFwdScore = FLT_MAX;
	m_Omega = FLT_MAX;
	m_OmegaFwd = FLT_MAX;
	m_MuPrefilterPatternStr = "?";
	m_MKFPatternStr = "?";
	}
