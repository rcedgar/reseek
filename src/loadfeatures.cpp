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
	asserta(optset_fdir && optset_fs);

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
		FEATURE F = LoadFeature(Path);
		double w = StrToFloat(Fields2[1]);
		AddFeature(F, w);
		ProgressLog("%s : %.3g\n", Path.c_str(), w);
		}
	SetScoreMxs();
	}
