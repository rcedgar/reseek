#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "sort.h"

vector<FEATURE> DSSParams::m_MuFeatures;
vector<uint> DSSParams::m_MuAlphaSizes;
uint DSSParams::m_MuAlphaSize = UINT_MAX;

void DSSParams::SetMuFeatures(const vector<FEATURE> &Fs)
	{
	m_MuAlphaSizes.clear();
	m_MuFeatures = Fs;
	m_MuAlphaSize = 1;
	for (uint i = 0; i < SIZE(Fs); ++i)
		{
		uint AS = DSS::GetAlphaSize(Fs[i]);
		m_MuAlphaSizes.push_back(AS);
		m_MuAlphaSize *= AS;
		}
	}

void DSSParams::SetFromCmdLine(uint DBSize)
	{
	int i = int(optset_fast) + int(optset_sensitive) + int(optset_verysensitive);
	if (i != 1)
		Die("Must specify -fast -sensitive or -verysensitive");

	if (optset_dbsize)
		m_DBSize = (float) opt_dbsize;
	else
		m_DBSize = (float) DBSize;

	m_Evalue_a = 4.0f;		if (optset_evalue_a) m_Evalue_a = float(opt_evalue_a);
	m_Evalue_b = -43.0f;	if (optset_evalue_b) m_Evalue_b = float(opt_evalue_b);

	vector<FEATURE> MuFeatures;
	MuFeatures.push_back(FEATURE_SS3);
	MuFeatures.push_back(FEATURE_NENSS3);
	MuFeatures.push_back(FEATURE_RENDist4);
	DSSParams::SetMuFeatures(MuFeatures);

	if (optset_namedparams)
		SetNamedParams(opt_namedparams);
	else if (optset_params)
		FromParamStr(opt_params);
	else if (optset_paramsf)
		FromTsv(opt_paramsf);
	else
		SetNamedParams("defaults");
	if (optset_fast)
		{
		m_Omega = 22;
		m_OmegaFwd = 50;
		m_MKFL = 500;
		m_MKF_X1 = 8;
		m_MKF_X2 = 8;
		m_MKF_MinHSPScore = 50;
		m_MKF_MinMegaHSPScore = -4;
		}
	else if (optset_sensitive)
		{
		m_Omega = 12;
		m_OmegaFwd = 20;
		m_MKFL = 600;
		m_MKF_X1 = 8;
		m_MKF_X2 = 8;
		m_MKF_MinHSPScore = 50;
		m_MKF_MinMegaHSPScore = -4;
		}
	else if (optset_verysensitive) 
		{
		m_Omega = 0;
		m_OmegaFwd = 0;
		m_MKFL = 99999;
		m_MKF_X1 = 99999;
		m_MKF_X2 = 99999;
		m_MKF_MinHSPScore = 0;
		m_MKF_MinMegaHSPScore = -99999;
		}
	const int MINUS = -1; // for visual emphasis here
	if (optset_omega) { m_Omega = (float) opt_omega;  }
	if (optset_omegafwd) { m_OmegaFwd = (float) opt_omegafwd; }
	if (optset_minfwdscore) { m_MinFwdScore = float(opt_minfwdscore); }
	if (optset_gapopen) { m_GapOpen =  MINUS*float(opt_gapopen); }
	if (optset_gapopen) { m_GapExt = MINUS*float(opt_gapext); }
	if (optset_para_mugapopen) { m_ParaMuGapOpen = opt_para_mugapopen; }
	if (optset_para_mugapext) { m_ParaMuGapExt = opt_para_mugapext; }
	if (optset_pattern) { m_PatternStr = string(opt_pattern); }
	if (optset_minhsp) { m_MKF_MinHSPScore = opt_minhsp; }
	if (optset_minmegahsp) { m_MKF_MinMegaHSPScore = float(opt_minmegahsp); }
	if (optset_xdrop1) { m_MKF_X1 = int(opt_xdrop1); }
	if (optset_xdrop2) { m_MKF_X2 = int(opt_xdrop2); }
	if (optset_mkfl) { m_MKFL = int(opt_mkfl); }

	if (m_GapOpen > 0 || m_GapExt > 0)
		Die("open=%.3g ext=%.3g, gap penalties must be >= 0",
		  opt_gapopen, opt_gapext);

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
		uint AS = g_AlphaSizes2[F];
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
	asserta(m_ScoreMxs == 0);
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

float DSSParams::GetEvalueGumbel(float TS, float mu, float beta) const
	{
	double gumbel_cdf(double mu, double beta, double x);
	double x = -log(TS);
	//double P = gumbel_cdf(2.5, 0.613, x);
	double P = gumbel_cdf(mu, beta, x);
	float Evalue = float(P*m_DBSize);
	return Evalue;
	}

float DSSParams::GetEvalueSlope(float TestStatistic, float m, float b) const
	{
	float PredMinusLogP = m*TestStatistic + b;
	float P = expf(-PredMinusLogP);
	float Evalue = P*m_DBSize;
	if (Evalue > 1)
		Evalue = log10f(Evalue) + 1;
	return Evalue;
	}

float DSSParams::GetEvalueOldLinear(float TestStatistic) const
	{
	const float Slope = m_Evalue_old_linear_Slope;
	const float Intercept = m_Evalue_linear_Intercept;
	float logNF = Slope*TestStatistic + Intercept;
	float NF = powf(10, logNF);
	float Evalue = NF*m_DBSize/1e8f;
	return Evalue;
	}

float DSSParams::GetEvalue(float TestStatistic) const
	{
	if (TestStatistic <= 0)
		return 99999;
	asserta(m_DBSize != 0 && m_DBSize != FLT_MAX);

	if (opt_gum)
		return GetEvalueGumbel(TestStatistic,
		  m_Evalue_Gumbel_mu, m_Evalue_Gumbel_beta);

	return GetEvalueSlope(TestStatistic,
	  m_Evalue_linear_m, m_Evalue_linear_b);
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
