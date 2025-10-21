#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "featuretrainer.h"

/***
Unweighted score mxs, overwrite for loaded features:
	trained_features.cpp

class DSSParams
	static bool m_ApplyWeightsDone;
	static vector<FEATURE> m_Features;
	static vector<float> m_Weights;
	static float ***m_ScoreMxs;			// weights built in

Converting float value to letter:
getbins.cpp(5)			void DSSParams::GetBins(), slow used only to initialize 
valuetoint_new.cpp(7)	static vector<vector<float> > s_BinTs(FEATURE_COUNT);
valuetoint_new.cpp(19):	static bool Init()
valuetoint_new.cpp(59)	uint DSSParams::ValueToInt_Feature(FEATURE F, float Value);
***/

//static vector<FEATURE> s_LoadedFeatures;
//static vector<float> s_LoadedWeights(FEATURE_COUNT);
//static vector<vector<float> > s_LoadedBinTs(FEATURE_COUNT);
//static vector<vector<vector<float> > > s_LoadedScoreMxs(FEATURE_COUNT);

FEATURE DSSParams::LoadFeature(const string &FN)
	{
	vector<vector<float> > ScoreMx;
	FILE *f = OpenStdioFile(FN);
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 2);
	FEATURE F = StrToFeature(Fields[0].c_str());
	uint AlphaSize = StrToUint(Fields[1]);
	ScoreMx.resize(AlphaSize);
	asserta(AlphaSize == GetAlphaSize(F));
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
	OverwriteUnweightedScoreMx(F, ScoreMx);

	if (!FeatureIsInt(F))
		{
		vector<float> BinTs;
		for (uint Letter = 0; Letter + 1 < AlphaSize; ++Letter)
			{
			bool Ok = ReadLineStdioFile(f, Line);
			asserta(Ok);
			Split(Line, Fields, '\t');
			asserta(SIZE(Fields) == 2);
			asserta(StrToUint(Fields[0]) == Letter);
			BinTs.push_back(StrToFloatf(Fields[1].c_str()));
			}
		OverwriteBinTs(F, BinTs);
		}
	CloseStdioFile(f);
	return F;
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
		asserta(SIZE(Fields) == 1);

		string Path = Fields[0];
		FEATURE F = LoadFeature(Path);
		ProgressLog("%s\n", Path.c_str());
		}
	ProgressLog("gapopen: %.3g\n", -m_GapOpen);
	ProgressLog("gapext: %.3g\n", -m_GapExt);
//	SetScoreMxs();
	}

void DSSParams::DumpFeature(FEATURE F, const string &FN)
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);
	uint AlphaSize = GetAlphaSize(F);
	fprintf(f, "%s\t%u\n", FeatureToStr(F), AlphaSize);
	for (uint Letter = 0; Letter < AlphaSize; ++Letter)
		{
		fprintf(f, "%u", Letter);
		for (uint Letter2 = 0; Letter2 < AlphaSize; ++Letter2)
			{
			float Score = g_ScoreMxs2[F][Letter][Letter2];
			fprintf(f, "\t%.4g", Score);
			}
		fprintf(f, "\n");
		}

	if (!FeatureIsInt(F))
		{
		vector<float> BinTs;
		GetBinTs(F, BinTs);
		asserta(SIZE(BinTs) + 1 == AlphaSize);
		for (uint i = 0; i + 1 < AlphaSize; ++i)
			fprintf(f, "%u\t%.3g\n", i, BinTs[i]);
		}
	CloseStdioFile(f);
	}

void cmd_load_features()
	{
	DSSParams::LoadFeatures(g_Arg1);
	DSSParams::SetFeatures();

	if (optset_output)
		{
		uint FeatureCount = DSSParams::GetFeatureCount();
		const string Prefix = opt(output);
		for (uint FIdx = 0; FIdx < FeatureCount; ++FIdx)
			{
			FEATURE F = DSSParams::m_Features[FIdx];
			const string &FN = Prefix + FeatureToStr(F);
			Progress("%s\n", FN.c_str());
			DSSParams::DumpFeature(F, FN);
			}
		}
	}
