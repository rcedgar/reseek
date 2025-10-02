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
