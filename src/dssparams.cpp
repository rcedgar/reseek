#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "sort.h"
#include "muprefilter_params.h"

const FEATURE DSSParams::m_MuFeatures[m_MuFeatureCount] =
	{
	FEATURE_SS3,
	FEATURE_NENSS3,
	FEATURE_RENDist4
	};
const uint DSSParams::m_MuAlphaSizes[m_MuFeatureCount] = {3, 3, 4};
uint const DSSParams::m_MuAlphaSize = 36;

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
	case DM_UseCommandLineOption:	return GetAlgoModeFromCommandLine(AM_Invalid);
	case DM_DefaultFast:			return GetAlgoModeFromCommandLine(AM_Fast);
	case DM_DefaultSensitive:		return GetAlgoModeFromCommandLine(AM_Sensitive);
		}
	asserta(false);
	return AM_Invalid;
	}

void DSSParams::SetDSSParams(DECIDE_MODE DM, uint DBSize)
	{
	SetDefaults();

	ALGO_MODE AM = GetAlgoMode(DM);

	if (DBSize == UINT_MAX)
		{
		if (optset_dbsize)
			m_DBSize = (float) opt(dbsize);
		else
			m_DBSize = SCOP40_DBSIZE;
		}
	else
		m_DBSize = (float) DBSize;

	switch (AM)
		{
	case AM_Fast:
		m_Omega = 22;
		m_OmegaFwd = 50;
		m_MKFL = 500;
		m_MKF_X1 = 8;
		//m_MKF_X2 = 8;
		m_MKF_MinHSPScore = 50;
		m_MKF_MinMegaHSPScore = -4;
		break;

	case AM_Sensitive:
		m_Omega = 12;
		m_OmegaFwd = 20;
		m_MKFL = 600;
		m_MKF_X1 = 8;
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
		break;

	default:
		asserta(false);
		}

	m_MKFPatternStr = "111";
	m_MuPrefilterPatternStr = string(muprefilter_pattern);

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
	}

void DSSParams::FromTsv(const string &FileName)
	{
	Clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	while (ReadLineStdioFile(f, Line))
		{
		vector<string> Fields;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);
		const string &Name = Fields[0];
		float Value = (float) StrToFloat(Fields[1]);
		SetParam(Name, Value, true);
		}
	CloseStdioFile(f);
	}

void DSSParams::ToFev(FILE *f, bool nl) const
	{
	if (f == 0)
		return;
	const uint FeatureCount = GetFeatureCount();
	fprintf(f, "NF=%u", FeatureCount);
	for (uint i = 0; i < FeatureCount; ++i)
		{
		FEATURE F = m_Features[i];
		fprintf(f, "\t%s=%.6g", FeatureToStr(F), m_Weights[i]);
		}
#define P(x)	fprintf(f, "\t%s=%.6g", #x, m_##x);
#include "scalarparams.h"

	if (nl)
		fprintf(f, "\n");
	}

uint DSSParams::GetFeatureCount() const
	{
	uint n = SIZE(m_Features);
	asserta(SIZE(m_Weights) == n);
	return n;
	}

float DSSParams::GetParam(const string &Name) const
	{
#define P(f)	if (Name == #f) { return m_##f; }
#include "scalarparams.h"

	for (uint F = 0; F < FEATURE_COUNT; ++F)
		{
		if (Name == FeatureToStr(F))
			{
			uint Idx = GetFeatureIdx(FEATURE(F));
			return m_Weights[Idx];
			}
		}
	Die("GetParam(%s)", Name.c_str());
	return FLT_MAX;
	}

int DSSParams::GetIntParam(const string &Name) const
	{
#define x(f)	if (Name == #f) { return m_##f; }
	x(ParaMuGapOpen);
	x(ParaMuGapExt);
#undef x
	Die("GetIntParam(%s)", Name.c_str());
	return INT_MAX;
	}

void DSSParams::SetIntParam(const string &Name, int Value)
	{
#define x(f)	if (Name == #f) { m_##f = Value; return; }
	x(ParaMuGapOpen);
	x(ParaMuGapExt);
#undef x
	Die("SetParam(%s)", Name.c_str());
	}

void DSSParams::SetParam(const string &Name, float Value, bool AppendIfWeight)
	{
#define P(f)	if (Name == #f) { m_##f = Value; return; }
#include "scalarparams.h"

	if (AppendIfWeight)
		{
		FEATURE F = StrToFeature(Name.c_str());
		m_Features.push_back(F);
		m_Weights.push_back(Value);
		return;
		}
	else
		{
		for (uint Idx = 0; Idx < SIZE(m_Features); ++Idx)
			{
			FEATURE F = m_Features[Idx];
			if (Name == FeatureToStr(F))
				{
				m_Weights[Idx] = Value;
				return;
				}
			}
		}
	Die("SetParam(%s)", Name.c_str());
	}

uint DSSParams::GetFeatureIdx(FEATURE F) const
	{
	for (uint Idx = 0; Idx < SIZE(m_Features); ++Idx)
		if (m_Features[Idx] == F)
			return Idx;
	Die("GetFeatureIdx(%u)", F);
	return UINT_MAX;
	}

uint DSSParams::GetFeatureIdx_NoError(FEATURE F) const
	{
	for (uint Idx = 0; Idx < SIZE(m_Features); ++Idx)
		if (m_Features[Idx] == F)
			return Idx;
	return UINT_MAX;
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

DSSParams::~DSSParams()
	{
	uint FeatureCount = GetFeatureCount();
	if (!m_OwnScoreMxs)
		return;
	for (uint Idx = 0; Idx < FeatureCount; ++Idx)
		{
		FEATURE F = m_Features[Idx];
		asserta(uint(F) < FEATURE_COUNT);
		uint AS = GetAlphaSize(F); // g_AlphaSizes2[F];
		for (uint Letter1 = 0; Letter1 < AS; ++Letter1)
			{
			asserta(m_ScoreMxs[F][Letter1] != 0);
			myfree(m_ScoreMxs[F][Letter1]);
			}
		asserta(m_ScoreMxs[F] != 0);
		myfree(m_ScoreMxs[F]);
		}
	myfree(m_ScoreMxs);
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
		uint AS = GetAlphaSize(F); // g_AlphaSizes2[F];
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
		uint AS = GetAlphaSize(F); // g_AlphaSizes2[F];
		if (AS == 0)
			Die("Feature %s not supported", FeatureToStr(F));
		m_ScoreMxs[F] = myalloc(float *, AS);
		for (uint Letter1 = 0; Letter1 < AS; ++Letter1)
			{
			const float * const *MxF = GetScoreMx(F);
			m_ScoreMxs[F][Letter1] = myalloc(float, AS);
			for (uint Letter2 = 0; Letter2 < AS; ++Letter2)
				m_ScoreMxs[F][Letter1][Letter2] = w*MxF[Letter1][Letter2]; // g_ScoreMxs2[F][Letter1][Letter2];
			}
		}
	}
