#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "sort.h"

double Q_func(double x, double mu, double sigma);

vector<FEATURE> DSSParams::m_ComboFeatures;
vector<uint> DSSParams::m_ComboAlphaSizes;
uint DSSParams::m_ComboAlphaSize = UINT_MAX;

void DSSParams::SetComboFeatures(const vector<FEATURE> &Fs)
	{
	m_ComboAlphaSizes.clear();
	m_ComboFeatures = Fs;
	m_ComboAlphaSize = 1;
	for (uint i = 0; i < SIZE(Fs); ++i)
		{
		uint AS = DSS::GetAlphaSize(Fs[i]);
		m_ComboAlphaSizes.push_back(AS);
		m_ComboAlphaSize *= AS;
		}
	}

void DSSParams::SetFromCmdLine(uint DBSize)
	{
	if (optset_dbsize)
		m_DBSize = (float) opt_dbsize;
	else
		m_DBSize = (float) DBSize;

	vector<FEATURE> ComboFeatures;
	ComboFeatures.push_back(FEATURE_SS3);
	ComboFeatures.push_back(FEATURE_NbrSS3);
	ComboFeatures.push_back(FEATURE_RevNbrDist4);
	DSSParams::SetComboFeatures(ComboFeatures);

	if (optset_namedparams)
		SetNamedParams(opt_namedparams);
	else if (optset_paramsf)
		FromTsv(opt_paramsf);
	else 
		SetNamedParams("defaults");

	const int MINUS = -1; // for visual emphasis here
	if (optset_omega) { m_Omega = (float) opt_omega; Psa(m_Desc, " -omega %.4g", opt_omega); }
	if (optset_daliw) { m_DALIw = (float) opt_daliw; Psa(m_Desc, " -daliw %.4g", opt_daliw); }
	if (optset_lambda) { m_Lambda = opt_lambda; Psa(m_Desc, " -lambda %u", opt_lambda); }
	if (optset_minfwdscore) { m_MinFwdScore = float(opt_minfwdscore); Psa(m_Desc, " -minfwdscore %.4g", opt_minfwdscore); }
	if (optset_gapopen) { m_GapOpen =  MINUS*float(opt_gapopen); Psa(m_Desc, " -gapopen %.4g", opt_gapopen); }
	if (optset_gapopen) { m_GapExt = MINUS*float(opt_gapext); Psa(m_Desc, " -gapext %.4g", opt_gapext); }
	if (optset_minu) { m_MinU = opt_minu; Psa(m_Desc, " -minu %u", opt_minu); }
	if (optset_maxaccepts) { m_MaxAccepts = opt_maxaccepts; Psa(m_Desc, " -maxaccepts %u", opt_maxaccepts); }
	if (optset_maxrejects) { m_MaxRejects = opt_maxrejects; Psa(m_Desc, " -maxrejects %u", opt_maxrejects); }
	if (optset_usort) { m_USort = true;  Psa(m_Desc, " -usort"); }
	if (optset_usecombopath) { m_UseComboPath = true; Psa(m_Desc, " -usecombopath"); Warning("-usecombopath bad idea"); }

	if (m_GapOpen > 0 || m_GapExt > 0)
		Die("open=%.3g ext=%.3g, gap penalties must be >= 0",
		  opt_gapopen, opt_gapext);

	InitScoreMxs();
	WriteSummary(g_fLog);
	WriteSummary(stderr);
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

void DSSParams::WriteSummary(FILE *f) const
	{
	if (f == 0)
		return;
	const uint FeatureCount = GetFeatureCount();
	fprintf(f, "---------------------------------------------------------------------------------\n");
	if (!m_Desc.empty())
		{
		fprintf(f, "%s\n", m_Desc.c_str());
		fprintf(f, "=================================================================================\n");
		}
	for (uint i = 0; i < FeatureCount; ++i)
		{
		FEATURE F = m_Features[i];
		if (i > 0 && i%5 == 0)
			fprintf(f, "\n");
		fprintf(f, "%s:%u/%.3f ",
		  FeatureToStr(F), DSS::GetAlphaSize(F), m_Weights[i]);
		}
	fprintf(f, "\n");
	fprintf(f, "GapO/E %.3f/", -m_GapOpen);
	fprintf(f, "%.3f", -m_GapExt);
	fprintf(f, " FwdM %.2f", m_FwdMatchScore);
	fprintf(f, " MinFS %.1f", m_MinFwdScore);
	fprintf(f, " Lamda %u", m_Lambda);
	if (m_Omega != FLT_MAX)
		fprintf(f, " Omega %.1f", m_Omega);
	fprintf(f, " DALIw %.1f", m_DALIw);
	fprintf(f, " MinU %u", m_MinU);
	fprintf(f, "\n");
	fprintf(f, "---------------------------------------------------------------------------------\n");
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
	x(ParaComboGapOpen);
	x(ParaComboGapExt);
#undef x
	Die("GetIntParam(%s)", Name.c_str());
	return INT_MAX;
	}

void DSSParams::SetIntParam(const string &Name, int Value)
	{
#define x(f)	if (Name == #f) { m_##f = Value; return; }
	x(ParaComboGapOpen);
	x(ParaComboGapExt);
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

// Superfamily: Linear fit to -log(P) m=20.5 b=2.89
// Fold:        Linear fit to -log(P) m=26.6 b=2.42
float DSSParams::GetEvalue(float TestStatistic, float m, float b) const
	{
	if (TestStatistic <= 0)
		return 99999;
	asserta(m_DBSize != 0 && m_DBSize != FLT_MAX);
	float Evalue = FLT_MAX;
	if (m != FLT_MAX)
		{
		asserta(b != FLT_MAX);
	//	PredMinusLogP = m*TS + b;
		float PredMinusLogP = m*TestStatistic + b;
		float P = expf(-PredMinusLogP);
		Evalue = P*m_DBSize;
		return Evalue;
		}

	if (opt_gum)
		{
		double gumbel_cdf(double mu, double beta, double x);
		double x = -log(TestStatistic);
		double P = gumbel_cdf(2.5, 0.613, x);
		Evalue = float(P*m_DBSize);
		}
	else
		{
		const float Slope = -6.6f;
		const float Intercept = 6.1f;
		float logNF = Slope*TestStatistic + Intercept;
		float NF = powf(10, logNF);
		float OldEvalue = NF*m_DBSize/1e8f;

		const float m = 20.5f;
		const float b = 2.90f;
		float MinusLogP = TestStatistic*m + b;
		float P = expf(-MinusLogP);
		Evalue = P*m_DBSize;
		//Log("%.3g\t%.3g\n", OldEvalue, Evalue);
		}
	return Evalue;
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
		m_ScoreMxs[F] = myalloc(float *, AS);
		for (uint Letter1 = 0; Letter1 < AS; ++Letter1)
			{
			m_ScoreMxs[F][Letter1] = myalloc(float, AS);
			for (uint Letter2 = 0; Letter2 < AS; ++Letter2)
				m_ScoreMxs[F][Letter1][Letter2] = w*g_ScoreMxs2[F][Letter1][Letter2];
			}
		}
	}
