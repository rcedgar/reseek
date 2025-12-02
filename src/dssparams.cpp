#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "sort.h"

vector<FEATURE> DSSParams::m_Features;
vector<float> DSSParams::m_Weights;

////////////////////////////////////////////////
// Mu-related parameters
////////////////////////////////////////////////
int DSSParams::m_ParaBits = 8;

int DSSParams::m_Omega8 = 22;
int DSSParams::m_OmegaFwd8 = 50;

int DSSParams::m_Omega16 = 250;
int DSSParams::m_OmegaFwd16 = 500;

int DSSParams::m_ParaMuGapOpen16 = 8;
int DSSParams::m_ParaMuGapExt16 = 5;
////////////////////////////////////////////////

float DSSParams::m_MinFwdScore = 7;
string DSSParams::m_MKFPatternStr =  "111";;
float ***DSSParams::m_ScoreMxs = 0;

int  DSSParams::m_PrefilterMinKmerPairScore = 36;
uint DSSParams::m_rsb_size = 1500;
uint DSSParams::m_MKFL = 500;
int DSSParams::m_MKF_X1 = 8;
int DSSParams::m_MKF_X2 = 8;
int DSSParams::m_MKF_MinMuHSPScore = 50;
float DSSParams::m_MKF_MinMegaHSPScore = -4;

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

float DSSParams::m_GapOpen = -0.685533f;
float DSSParams::m_GapExt = -0.051881f;

float DSSParams::m_dpw = 1.7f;
float DSSParams::m_lddtw = 0.13f;
float DSSParams::m_ladd = 250.0f;
float DSSParams::m_revtsw = 2.0f;

// Used to initialize g_AlphaSizes2, careful of
//	order dependencies in static bool Init() idiom.
uint DSSParams::GetAlphaSize(FEATURE F, bool FailOk)
	{
	switch (F)
		{
	case FEATURE_AA:
	case FEATURE_B62:
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
	case FEATURE_PENDist:
	case FEATURE_MENDist:
	case FEATURE_PENConf:
	case FEATURE_MENConf:
	case FEATURE_PMDist:
	case FEATURE_PMDistDiff:
	case FEATURE_HelixDens:
	case FEATURE_StrandDens:
	case FEATURE_DstNxtHlx:
	case FEATURE_DstPrvHlx:
	case FEATURE_NX:
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

int DSSParams::GetOmega()
	{
	if (m_ParaBits == 8)
		return m_Omega8;
	else if (m_ParaBits == 16)
		return m_Omega16;
	Die("GetOmega");
	return false;
	}

int DSSParams::GetOmegaFwd()
	{
	if (m_ParaBits == 8)
		return m_OmegaFwd8;
	else if (m_ParaBits == 16)
		return m_OmegaFwd16;
	Die("GetOmega");
	return false;
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
	X(Omega8,				22,		12,		0)
	X(OmegaFwd8,			50,		20,		0)
	X(Omega16,				16,		12,		0)
	X(OmegaFwd16,			48,		20,		0)
	X(MKFL,					500,	600,	99999)
	X(MKF_X1,				8,		8,		99999)
	X(MKF_X2,				8,		8,		99999)
	X(MKF_MinMuHSPScore,	50,		50,		0)
	X(MKF_MinMegaHSPScore,	-4,		-4,		-99999)
#undef X

	if (optset_rsb_size)	m_rsb_size = opt(rsb_size);
	if (optset_parabits)	m_ParaBits = opt(parabits);
	if (optset_omega8)		m_Omega8 = (int) opt(omega8);
	if (optset_omega16)		m_Omega16 = (int) opt(omega16);
	if (optset_omegafwd8)	m_OmegaFwd8 = (int) opt(omegafwd8);
	if (optset_omegafwd16)	m_OmegaFwd16 = (int) opt(omegafwd16);
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
	AllocScoreMxs();
	NormalizeWeights();
	CreateWeightedScoreMxs();
	}

void DSSParams::FreeScoreMxs()
	{
	if (m_ScoreMxs == 0)
		return;

	for (uint i = 0; i < SIZE(m_Features); ++i)
		{
		FEATURE F = m_Features[i];
		if (m_ScoreMxs[F] != 0)
			{
			uint AS = GetAlphaSize(F);
			for (uint j = 0; j < AS; ++j)
				{
				asserta(m_ScoreMxs[F][j] != 0);
				myfree(m_ScoreMxs[F][j]);
				}
			myfree(m_ScoreMxs[F]);
			m_ScoreMxs[F] = 0;
			}
		}

	//for (uint i = 0; i < FEATURE_COUNT; ++i)
	//	asserta(m_ScoreMxs[i] == 0);

	myfree(m_ScoreMxs);
	}

void DSSParams::AllocScoreMxs()
	{
	if (m_ScoreMxs != 0)
		FreeScoreMxs();
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

void DSSParams::SetDefaultNonFeatureTunableParams()
	{
	m_GapOpen = -0.685533f;
	m_GapExt = -0.051881f;

	m_dpw = 1.7f;
	m_lddtw = 0.13f;
	m_ladd = 250.0f;
	m_revtsw = 2.0f;
	}


// Assume defaults are already set, except for features+weights
// Set features+weights from Str
// Overwrite any others in Str
void DSSParams::SetParamsFromStr(const string &Str)
	{
	m_Features.clear();
	m_Weights.clear();

	vector<string> Fields;
	Split(Str, Fields, ';');

	const uint n = SIZE(Fields);
	for (uint i = 0; i < n; ++i)
		{
		const string &NameEqValue = Fields[i];
		vector<string> Fields2;
		Split(NameEqValue, Fields2, '=');
		if (SIZE(Fields2) != 2)
			Die("DSSParams::SetParamsFromStr(%s) not name=value '%s'",
				Str.c_str(), Fields[i].c_str());
		const string &Name = Fields2[0];
		const string &ValueStr = Fields2[1];
		FEATURE F = StrToFeature(Name.c_str(), true);
		if (F == FEATURE(UINT_MAX))
			SetParam(Name, ValueStr);
		else
			{
			float Weight = StrToFloatf(ValueStr);
			m_Features.push_back(F);
			m_Weights.push_back(Weight);
			}
		}
	asserta(SIZE(m_Features) > 0);
	SetScoreMxsFromFeatures();
	}

void DSSParams::OverwriteFeatures(const vector<FEATURE> &Fs,
	const vector<float> &Weights)
	{
	m_Features.clear();
	m_Weights.clear();
	const uint n = SIZE(Fs);
	asserta(SIZE(Weights) == n);
	for (uint i = 0; i < n; ++i)
		{
		m_Features.push_back(Fs[i]);
		m_Weights.push_back(Weights[i]);
		}
	SetScoreMxsFromFeatures();
	}

void DSSParams::SetStandardFeatures()
	{
	if (opt(newparams))
		{
		// git [3d7ef20]
		// scop40c SEPQ0.1=0.282 SEPQ1=0.374 SEPQ10=0.452 Area0=0.979 Sum3=1.577 -fast ../data/scop40c.bca
		// scop40  SEPQ0.1=0.226 SEPQ1=0.330 SEPQ10=0.423 Area0=0.829 Area3=1.156 Sum3=1.370
		// AA=4.4E-01;NENDist=1.6E-01;Conf=1.9E-01;NENConf=8.3E-02;RENDist=5.2E-02;DstNxtHlx=4.4E-02;StrandDens=0.0E+00;NormDens=4.0E-02;
		AddFeature(FEATURE_AA,			4.4E-01f);
		AddFeature(FEATURE_NENDist,		1.6E-01f);
		AddFeature(FEATURE_Conf,		1.9E-01f);
		AddFeature(FEATURE_NENConf,		8.3E-02f);
		AddFeature(FEATURE_RENDist,		5.2E-02f);
		AddFeature(FEATURE_DstNxtHlx,	4.4E-02f);
		AddFeature(FEATURE_NormDens,	4.0E-02f);

		// gap2=7.6E-01;dpw=2.0E+00;lddtw=2.0E-01;revtsw=2.5E+00;logladd=2.4E+00;
		m_GapOpen = -7.67E-01f;
		m_GapExt = m_GapOpen/10;
		m_dpw = 2.0E+00f;
		m_lddtw = 2.0E-01f;
		m_revtsw = 2.5E+00f;
		m_ladd = powf(10, 2.4E+00f);
		}
	else
		{
		// scop40c SEPQ0.1=0.276 SEPQ1=0.368 SEPQ10=0.449 Area0=0.960             Sum3=1.553 -fast ../data/scop40c.bca
		// scop40  SEPQ0.1=0.211 SEPQ1=0.324 SEPQ10=0.423 Area0=0.800 Area3=1.119 Sum3=1.351
		AddFeature(FEATURE_AA,			0.398145f);
		AddFeature(FEATURE_NENDist,		0.129367f);
		AddFeature(FEATURE_Conf,		0.202354f);
		AddFeature(FEATURE_NENConf,		0.149383f);
		AddFeature(FEATURE_RENDist,		0.0937677f);
		AddFeature(FEATURE_DstNxtHlx,	0.00475462f);
		AddFeature(FEATURE_StrandDens,	0.0183853f);
		AddFeature(FEATURE_NormDens,	0.00384384f);
		}
	SetScoreMxsFromFeatures();
	}

void DSSParams::Init(DECIDE_MODE DM)
	{
	if (optset_feature_spec)
		LoadFeatures(opt(feature_spec));
	else if (optset_params)
		SetParamsFromStr(opt(params));
	else
		SetStandardFeatures();
	SetAlgoMode(DM);
	string ParamsStr;
	DSSParams::GetParamsStr(ParamsStr);
	Log("DSSParams::Init(%d) %s\n", int(DM), ParamsStr.c_str());
	}

void DSSParams::LogMe()
	{
	Log("%u features\n", SIZE(m_Features));
	for (uint i = 0; i < SIZE(m_Features); ++i)
		Log("%s %.3g\n", FeatureToStr(m_Features[i]), m_Weights[i]);

	LogParamData();
	Log("m_MKFPatternStr %s\n", m_MKFPatternStr.c_str());
	}

const float * const *DSSParams::GetScoreMx(FEATURE F)
	{
	return g_ScoreMxs2[F];
	}
