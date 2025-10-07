#include "myutils.h"
#include "dssparams.h"
#include "featuretrainer.h"
#include "dss.h"

extern float **g_ScoreMxs2[FEATURE_COUNT];
extern float **g_FreqMxs2[FEATURE_COUNT];
extern float *g_FreqVecs2[FEATURE_COUNT];
extern uint g_AlphaSizes2[FEATURE_COUNT];

static void DumpOldFeature(const string &OutPrefix, FEATURE F)
	{
	const string &FeatureName = string(FeatureToStr(F));
	const string FN = OutPrefix + FeatureName + ".tsv";
	FILE *f = CreateStdioFile(FN);
	uint AS = g_AlphaSizes2[F];
	asserta(AS > 0);
	bool IsInt = FeatureIsInt(F);

//feature NENDist
//type    float
//undef_binning   zeroov
//min     2.5
//med     6.86
//max     48
//undef   0.0000
//alpha_size      16
//expected_score  0.838

	fprintf(f, "feature\t%s\n", FeatureName.c_str());
	fprintf(f, "type\t%s\n", IsInt ? "int" : "float");
	if (IsInt)
		fprintf(f, "undef_binning\tint\n");
	else
		{
		fprintf(f, "undef_binning\tzeroov\n");
		fprintf(f, "min\t0\n");
		fprintf(f, "med\t999\n");
		fprintf(f, "max\t9999\n");
		fprintf(f, "undef\t0.1234\n");
		}
	fprintf(f, "alpha_size\t%u\n", AS);
	fprintf(f, "expected_score\t0.555\n");

	const float *Freqs = g_FreqVecs2[F];
	float **FreqMx = g_FreqMxs2[F];

	const uint N = 1024*1024*1024;
	vector<uint> Counts;
	vector<vector<uint> > CountMx(AS);
	for (uint Letter1 = 0; Letter1 < AS; ++Letter1)
		{
		CountMx[Letter1].resize(AS);
		uint LetterCount = 0;
		for (uint Letter2 = 0; Letter2 < AS; ++Letter2)
			{
			uint n = uint(FreqMx[Letter1][Letter2]*N + 0.5);
			CountMx[Letter1][Letter2] = n;
			LetterCount += n;
			}
		Counts.push_back(LetterCount);
		}

	for (uint Letter = 0; Letter < AS; ++Letter)
		fprintf(f, "count\t%u\t%u\n", Letter, Counts[Letter]);

	for (uint Letter1 = 0; Letter1 < AS; ++Letter1)
		{
		fprintf(f, "pair_counts\t%u", Letter1);
		for (uint Letter2 = 0; Letter2 < AS; ++Letter2)
			{
			uint Count = CountMx[Letter1][Letter2];
			fprintf(f, "\t%u", Count);
			}
		fprintf(f, "\n");
		}
	if (!IsInt)
		{
		vector<float> Bins;
		DSS::GetBins(F, Bins);
		for (uint i = 0; i < SIZE(Bins); ++i)
			fprintf(f, "bint\t%u\t%.4g\n", i, Bins[i]);
		}

#if 1
//@@TODO dumps old ScoreMx, not consistent with counts
	float **OldScoreMx = g_ScoreMxs2[F];
	for (uint Letter1 = 0; Letter1 < AS; ++Letter1)
		{
		fprintf(f, "scoremx\t%u", Letter1);
		for (uint Letter2 = 0; Letter2 < AS; ++Letter2)
			fprintf(f, "\t%.5g", OldScoreMx[Letter1][Letter2]);
		fprintf(f, "\n");
		}

#else
	FeatureTrainer FT;
	vector<float> Freqs2;
	vector<vector<float> > FreqMx2;
	vector<vector<float> > ScoreMx2;

	const string FN2 = OutPrefix + FeatureName + ".scoremx";

	FILE *f2 = CreateStdioFile(FN2);
	FT.FromTsv(FN);
	FT.GetFreqs(Freqs2);
	FT.GetFreqMx(ScoreMx2);
	FT.GetLogOddsMx(ScoreMx2);
	FT.ScoreMxToSrc(f2);
	CloseStdioFile(f2);
#endif 
	
	CloseStdioFile(f);
	}

void cmd_dump_features()
	{
	DSSParams::SetDefaults();
	uint FeatureCount = DSSParams::GetFeatureCount();

	string OutPrefix = g_Arg1;
	for (uint i = 0; i < FeatureCount; ++i)
		{
		FEATURE F = DSSParams::m_Features[i];
		DumpOldFeature(OutPrefix, F);
		}
	}
