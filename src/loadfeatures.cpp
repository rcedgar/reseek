#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "featuretrainer.h"

FEATURE DSSParams::LoadNewTrainFeature(const string &FN)
	{
	FeatureTrainer FT;
	vector<float> Freqs;
	vector<vector<float> > FreqMx;
	vector<vector<float> > ScoreMx;
	FT.FromTsv(FN);
	FT.GetFreqs(Freqs);
	FT.GetFreqMx(FreqMx);
	FT.GetLogOddsMx(ScoreMx);
	DSS::SetNewTrainFeature(FT.m_F, Freqs, FreqMx, ScoreMx, FT.m_BinTs);
	return FT.m_F;
	}

void DSSParams::LoadNewTrainFeaturesFromCmdLine()
	{
	Clear();
	SetDefaults_Other();

	string FDir = opt(fdir);
	if (FDir == "")
		FDir = ".";
	Dirize(FDir);

	const string &fs = opt(fs);
	vector<string> Fields;
	Split(fs, Fields, '_');

	vector<string> FNs;
	vector<float> Weights;
	const uint N = SIZE(Fields);
	for (uint i = 0; i < N; ++i)
		{
		const string &Field = Fields[i];
		vector<string> Fields2;
		Split(Field, Fields2, ':');
		asserta(SIZE(Fields2) == 2);
		const string &FN = Fields2[0];
		string Path = FDir + FN;
		FEATURE F = LoadNewTrainFeature(Path);
		double w = StrToFloat(Fields2[1]);
		AddFeature(F, w);
		ProgressLog("%s : %.3g\n", Path.c_str(), w);
		}
	m_NewTrain = true;
	}
