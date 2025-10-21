#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "featuretrainer.h"

static vector<FEATURE> s_LoadedFeatures;
static vector<float> s_LoadedWeights(FEATURE_COUNT);
static vector<vector<float> > s_LoadedBinTs(FEATURE_COUNT);
static vector<vector<vector<float> > > s_LoadedScoreMxs(FEATURE_COUNT);

void DSSParams::LoadFeature(const string &FN,
	FEATURE &F, uint &AlphaSize, float &Weight,
	vector<vector<float> > &ScoreMx, vector<float> &BinTs)
	{
	ScoreMx.clear();
	BinTs.clear();
	ScoreMx.resize(AlphaSize);
	FILE *f = OpenStdioFile(FN);
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 3);
	F = StrToFeature(Fields[0].c_str());
	AlphaSize = StrToUint(Fields[1]);
	Weight = StrToFloatf(Fields[2]);

	for (uint Letter = 0; Letter < AlphaSize; ++Letter)
		{
		ScoreMx[Letter].resize(AlphaSize, FLT_MAX);
		bool Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == AlphaSize+1);
		asserta(StrToUint(Fields[0]) == Letter);
		for (uint Letter2 = 0; Letter2 < AlphaSize; ++Letter2)
			ScoreMx[Letter][Letter2] = StrToFloatf(Fields[Letter2+1]);
		}

	if (!FeatureIsInt(F))
		{
		for (uint Letter = 0; Letter < AlphaSize; ++Letter)
			{
			bool Ok = ReadLineStdioFile(f, Line);
			asserta(Ok);
			Split(Line, Fields, '\t');
			asserta(SIZE(Fields) == 2);
			asserta(StrToUint(Fields[0]) == Letter);
			BinTs.push_back(StrToFloatf(Fields[1].c_str()));
			}
		}

	CloseStdioFile(f);
	}

void DSSParams::LoadFeatures(const string &FN)
	{
	vector<string> Lines;
	vector<string> Fields;
	ReadLinesFromFile(FN, Lines);
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
		float w = StrToFloatf(Fields[1]);

		FEATURE F;
		vector<vector<float> > ScoreMx;
		vector<float> BinTs;
		uint AlphaSize;
		float Weight;
		LoadFeature(Path, F, AlphaSize, Weight, ScoreMx, BinTs);

		s_LoadedFeatures.push_back(F);
		s_LoadedWeights[F] = Weight;
		s_LoadedBinTs[F] = BinTs;
		s_LoadedScoreMxs[F] = ScoreMx;

		AddFeature(F, w);
		ProgressLog("%s : %.3g\n", Path.c_str(), w);
		}
	ProgressLog("gapopen: %.3g\n", -m_GapOpen);
	ProgressLog("gapext: %.3g\n", -m_GapExt);
//	SetScoreMxs();
	}

void cmd_load_features()
	{
	DSSParams::LoadFeatures(g_Arg1);
	uint FeatureCount = DSSParams::GetFeatureCount();

	//string OutPrefix = g_Arg1;
	//for (uint i = 0; i < FeatureCount; ++i)
	//	{
	//	FEATURE F = Params.m_Features[i];
	//	DumpFeature(OutPrefix, F);
	//	}
	}
