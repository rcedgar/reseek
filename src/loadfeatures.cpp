#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "featuretrainer.h"

FEATURE DSSParams::LoadFeature(const string &FN)
	{
	FeatureTrainer FT;
	vector<float> Freqs;
	vector<vector<float> > FreqMx;
	vector<vector<float> > ScoreMx;
	FT.FromTsv(FN);
	FT.GetFreqs(Freqs);
	FT.GetFreqMx(FreqMx);
	FT.GetLogOddsMx(ScoreMx);
	DSS::SetFeature(FT.m_F, FT.m_UB, Freqs, FreqMx,
					ScoreMx, FT.m_BinTs);
	return FT.m_F;
	}

void DSSParams::LoadFeatures()
	{
	asserta(optset_feature_spec);

	Clear();
	SetDefaults_Other();

	vector<string> Lines;
	vector<string> Fields;
	ReadLinesFromFile(opt(feature_spec), Lines);
	const uint N = SIZE(Lines);

	vector<string> FNs;
	vector<float> Weights;
	for (uint i = 0; i < N; ++i)
		{
		string &Line = Lines[i];
		StripWhiteSpace(Line);
		if (Line == "")
			continue;
		if (StartsWith(Line, "#"))
			continue;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);

		string Path = Fields[0];
		double w = StrToFloat(Fields[1]);

		FEATURE F = LoadFeature(Path);
		AddFeature(F, w);
		ProgressLog("%s : %.3g\n", Path.c_str(), w);
		}
	ProgressLog("gapopen: %.3g\n", -m_GapOpen);
	ProgressLog("gapext: %.3g\n", -m_GapExt);
	SetScoreMxs();
	}

void DSSParams::SetDefaults_Other()
	{
	m_GapOpen = -0.685533f;
	m_GapExt = -0.051881f;
	m_MinFwdScore = 7.0f;
	m_MuPrefilterPatternStr = "1110011";
	m_MKFPatternStr = "111";
	}

//void DSSParams::SetDefaults()
//	{
//	Clear();
//	SetDefaults_Features();
//	SetDefaults_Other();
//	}

void DSSParams::SetDefaults_Features()
	{
	Die("Default features not supported");

	m_Features.clear();
	m_Weights.clear();

	AddFeature(FEATURE_AA,			0.398145);
	AddFeature(FEATURE_NENDist,		0.129367);
	AddFeature(FEATURE_Conf,		0.202354);
	AddFeature(FEATURE_NENConf,		0.149383);
	AddFeature(FEATURE_RENDist,		0.0937677);
	AddFeature(FEATURE_DstNxtHlx,	0.00475462);
	AddFeature(FEATURE_StrandDens,	0.0183853);
	AddFeature(FEATURE_NormDens,	0.00384384);
	SetScoreMxs();
	}
