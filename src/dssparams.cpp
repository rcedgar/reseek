#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "sort.h"

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

int DSSParams::m_Density_W = 50;
int DSSParams::m_Density_w = 3;
int DSSParams::m_SSDensity_W = 50;
int DSSParams::m_SSDensity_w = 8;
float DSSParams::m_Density_Radius = 20.0;
float DSSParams::m_NU_ND_Radius = 20.0;
int DSSParams::m_NEN_W = 100;
int DSSParams::m_NEN_w = 12;
int DSSParams::m_NUDX_W = 50;
float DSSParams::m_DefaultNENDist = 10.0;
float DSSParams::m_SSDensity_epsilon = 1;
uint DSSParams::m_SSE_MinLength = 8;
uint DSSParams::m_SSE_Margin = 8;
uint DSSParams::m_PMDelta = 8;

// Used to initialize g_AlphaSizes2, careful of
//	order dependencies in static bool Init() idiom.
uint DSSParams::GetAlphaSize(FEATURE F, bool FailOk)
	{
	switch (F)
		{
	case FEATURE_AA:
		return 20;

	case FEATURE_SS:
	case FEATURE_NENSS:
	case FEATURE_RENSS:
	case FEATURE_NormDens4:
	case FEATURE_NENDist4:
	case FEATURE_RENDist4:
	case FEATURE_AA4:
		return 4;

	case FEATURE_SS3:
	case FEATURE_NENSS3:
	case FEATURE_RENSS3:
	case FEATURE_AA3:
		return 3;

	case FEATURE_Conf:
	case FEATURE_NENConf:
	case FEATURE_RENConf:
	case FEATURE_NormDens:
	case FEATURE_NENDist:
	case FEATURE_RENDist:
	case FEATURE_HelixDens:
	case FEATURE_StrandDens:
	case FEATURE_DstNxtHlx:
	case FEATURE_DstPrvHlx:
	case FEATURE_NX:
	case FEATURE_PMDist:
		return 16;

	case FEATURE_Mu:
		return 36;

	case FEATURE_ConfU:
		return 17;
		}
	if (!FailOk)
		Die("DSSParams::GetAlphaSize(%d)", int(F));
	return UINT_MAX;
	}

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
	asserta(m_ScoreMxs == 0);
	AllocScoreMxs();
	if (optset_params)
		SetParamsFromStr(opt(params));
	NormalizeWeights();
	CreateWeightedScoreMxs();
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
		uint AS = DSSParams::GetAlphaSize(F);
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
	}

void DSSParams::UpdateWeights(const vector<float> &Weights)
	{
	asserta(m_ScoreMxs != 0);
	asserta(SIZE(Weights) == SIZE(m_Weights));
	m_Weights = Weights;
	NormalizeWeights();
	CreateWeightedScoreMxs();
	}

void DSSParams::CreateWeightedScoreMxs()
	{
	asserta(m_ScoreMxs != 0);
	uint FeatureCount = GetFeatureCount();
	asserta(SIZE(m_Features) == FeatureCount);
	asserta(SIZE(m_Weights) == FeatureCount);
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
				{
				float Score = g_ScoreMxs2[F][Letter1][Letter2];
				m_ScoreMxs[F][Letter1][Letter2] = w*Score;
				}
			}
		}
	}

void DSSParams::SetParamStr(const string &Name, const string &StrValue)
	{
	FEATURE F = StrToFeature(Name.c_str(), true);
	float Value = StrToFloatf(StrValue);
	if (F != UINT_MAX)
		{
		const uint FeatureCount = GetFeatureCount();
		asserta(SIZE(m_Weights) == FeatureCount);
		for (uint k = 0; k < FeatureCount; ++k)
			{
			if (m_Features[k] == F)
				{
				m_Weights[k] = Value;
				return;
				}
			}
		Die("SetParamStr(%s, %s) not active feature",
			Name.c_str(), StrValue.c_str());
		}
	if (Name == "open") { m_GapOpen = -Value; return; }
	if (Name == "ext") { m_GapExt = -Value; return; }
	if (Name == "gap") { m_GapOpen = m_GapExt = -Value; return; }
	if (Name == "gap2") { m_GapOpen = -Value; m_GapExt = -Value*0.1f; }

	Die("SetParamStr(%s, %s) unknown param",
		Name.c_str(), StrValue.c_str());
	}

void DSSParams::GetParamStr(string &Str)
	{
	Str.clear();
	const uint FeatureCount = GetFeatureCount();
	asserta(SIZE(m_Features) == FeatureCount);
	asserta(SIZE(m_Weights) == FeatureCount);
	for (uint k = 0; k < FeatureCount; ++k)
		Psa(Str, "%s=%.3g;",
			FeatureToStr(m_Features[k]),
			m_Weights[k]);
	Psa(Str, "open=%.3g;", -m_GapOpen);
	Psa(Str, "ext=%.3g;", -m_GapExt);
	}

void DSSParams::SetParamsFromStr(const string &Str)
	{
	vector<string> Fields;
	Split(Str, Fields, ';');
	const uint n = SIZE(Fields);
	for (uint i = 0; i < n; ++i)
		{
		const string &NameEqValue = Fields[i];
		vector<string> Fields2;
		Split(NameEqValue, Fields2, '=');
		if (SIZE(Fields2) != 2)
			Die("Expected name=value in field %u '%s'", i, Str.c_str());
		const string &Name = Fields2[0];
		const string &Value = Fields2[1];
		SetParamStr(Name, Value);
		}
	}

void DSSParams::SetStandardFeatures()
	{
	AddFeature(FEATURE_AA,			0.398145f);
	AddFeature(FEATURE_NENDist,		0.129367f);
	AddFeature(FEATURE_Conf,		0.202354f);
	AddFeature(FEATURE_NENConf,		0.149383f);
	AddFeature(FEATURE_RENDist,		0.0937677f);
	AddFeature(FEATURE_DstNxtHlx,	0.00475462f);
	AddFeature(FEATURE_StrandDens,	0.0183853f);
	AddFeature(FEATURE_NormDens,	0.00384384f);
	if (optset_params)
		SetParamsFromStr(opt(params));
	SetScoreMxsFromFeatures();
	}

void DSSParams::Init(DECIDE_MODE DM)
	{
	if (optset_feature_spec)
		LoadFeatures(opt(feature_spec));
	else
		SetStandardFeatures();
	SetAlgoMode(DM);
	if (optset_gapopen) DSSParams::m_GapOpen = (float) opt(gapopen);
	if (optset_gapext) DSSParams::m_GapExt = (float) opt(gapext);
	if (optset_gap2) { DSSParams::m_GapOpen = DSSParams::m_GapExt = (float) opt(gap2); }
	string ParamStr;
	DSSParams::GetParamStr(ParamStr);
	Log("DSSParams::Init(%d) %s\n", int(DM), ParamStr.c_str());
	}