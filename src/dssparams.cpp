#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "sort.h"
#include "prefiltermuparams.h"

bool DSSParams::m_InitDone = false;
bool DSSParams::m_ApplyWeightsDone = false;
vector<FEATURE> DSSParams::m_Features;
vector<float> DSSParams::m_Weights;

float DSSParams::m_MinFwdScore = FLT_MAX;
float DSSParams::m_Omega = FLT_MAX;
float DSSParams::m_OmegaFwd = FLT_MAX;
string DSSParams::m_MKFPatternStr =  "111";;
float ***DSSParams::m_ScoreMxs = 0;

int  DSSParams::m_PrefilterMinKmerPairScore = 36;
uint DSSParams::m_rsb_size = 1500;
uint DSSParams::m_MKFL = UINT_MAX;
int DSSParams::m_MKF_X1 = INT_MAX;
int DSSParams::m_MKF_X2 = INT_MAX;
int DSSParams::m_MKF_MinMuHSPScore = INT_MAX;
float DSSParams::m_MKF_MinMegaHSPScore = FLT_MAX;

float DSSParams::m_GapOpen = -0.685533f;
float DSSParams::m_GapExt = -0.051881f;
int DSSParams::m_ParaMuGapOpen = 2;
int DSSParams::m_ParaMuGapExt = 1;

float DSSParams::m_dpw = 1.7f;
float DSSParams::m_lddtw = 0.13f;
float DSSParams::m_ladd = 250.0f;
float DSSParams::m_revtsw = 2.0f;

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

void DSSParams::SetAlgoMode(DECIDE_MODE DM)
	{
	ALGO_MODE AM = GetAlgoMode(DM);

#define X(x, xfast, xsens, xvsens)  \
	switch (AM)  \
		{  \
		case AM_Fast:			m_##x = xfast; break;  \
		case AM_Sensitive:		m_##x = xsens; break;  \
		case AM_VerySensitive:	m_##x = xvsens; break;  \
		default: asserta(false); \
		}

	X(MinFwdScore,			7.0f,	 7.0f,	0)
	X(Omega,				22,		12,		0)
	X(OmegaFwd,				50,		20,		0)
	X(MKFL,					500,	600,	99999)
	X(MKF_X1,				8,		8,		99999)
	X(MKF_X2,				8,		8,		99999)
	X(MKF_MinMuHSPScore,	50,		50,		0)
	X(MKF_MinMegaHSPScore,	-4,		-4,		-99999)
#undef X

	if (optset_rsb_size)
		m_rsb_size = opt(rsb_size);
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

void DSSParams::SetScoreMxsFromFeatures()
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
	}

void DSSParams::ApplyWeights()
	{
	asserta(m_ScoreMxs != 0);
	asserta(!m_ApplyWeightsDone);
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
	m_ApplyWeightsDone = true;
	}

void DSSParams::SetFeatures()
	{
	AddFeature(FEATURE_AA,			0.398145);
	AddFeature(FEATURE_NENDist,		0.129367);
	AddFeature(FEATURE_Conf,		0.202354);
	AddFeature(FEATURE_NENConf,		0.149383);
	AddFeature(FEATURE_RENDist,	0.0937677);
	AddFeature(FEATURE_DstNxtHlx,	0.00475462);
	AddFeature(FEATURE_StrandDens,	0.0183853);
	AddFeature(FEATURE_NormDens,	0.00384384);
	SetScoreMxsFromFeatures();
	}

void DSSParams::Init(DECIDE_MODE DM)
	{
	SetFeatures();
	SetAlgoMode(DM);
	}