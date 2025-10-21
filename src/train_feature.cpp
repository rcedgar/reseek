#include "myutils.h"
#include "featuretrainer.h"
#include "sort.h"

void cmd_train_feature()
	{
	asserta(optset_feature);
	DSSParams::Init(DM_DefaultFast);

	opt_force_undef = true;
	optset_force_undef = true;

	const string FeatureName = string(opt(feature));
	FEATURE F = StrToFeature(FeatureName.c_str());

	DSSParams::Init(DM_DefaultSensitive);

	uint AlphaSize = UINT_MAX;
	bool IsInt = FeatureIsInt(F);
	if (IsInt)
		AlphaSize = DSSParams::GetAlphaSize(F);
	else
		{
		asserta(optset_alpha_size);
		AlphaSize = opt(alpha_size);
		}

	FeatureTrainer FT;
	FT.ReadChains(g_Arg1);
	FT.ReadAlns(opt(input), true);
	if (optset_input2)
		FT.ReadAlns(opt(input2), false);
	ProgressLog("%u TPs, %u FPs\n", FT.m_TPCount, FT.m_FPCount);
	FT.SetFeature(F, AlphaSize);
	FT.Train(string(opt(style)));
	FT.SetAlnSubstScores();
	//return;//@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	float MinGapOpen = -0.2f;
	float MinGapExt = -0.2f;
	float dGapOpen = 0.2f;
	float dGapExt = 0.05f;
	float SimpleGridBestArea = 0;
	for (uint i = 0; i < 10; ++i)
		for (uint j = 0; j < 10; ++j)
			{
			float GapOpen = MinGapOpen + i*dGapOpen;
			float GapExt = MinGapExt + j*dGapExt;
			FT.SetAlnScoresAndArea(GapOpen, GapExt);
			SimpleGridBestArea = max(FT.m_Area, SimpleGridBestArea);
			}
	FT.SetAlnScoresAndArea(0, 0);//@@@@@@@@@@
	FT.LogROCStepsAndArea();

	vector<float> Areas;
	Areas.push_back(SimpleGridBestArea);
	const int ITERS = 8;
	for (int Iter = 0; Iter < ITERS; ++Iter)
		{
		FT.OptimizeGapPenalties();
		Areas.push_back(FT.m_Area);
		}

	vector<uint> Order(SIZE(Areas));
	QuickSortOrderDesc(Areas.data(), SIZE(Areas), Order.data());
	ProgressLog("Areas: ");
	for (uint Iter = 0; Iter < SIZE(Areas); ++Iter)
		{
		uint k = Order[Iter];
		ProgressLog(" %.3f", Areas[k]);
		if (Areas[k] == SimpleGridBestArea)
			ProgressLog("*");
		}
	ProgressLog("\n");

	FT.ToTsv(opt(output));
	vector<vector<float> > ScoreMx;
	float ES = FT.GetLogOddsMx(ScoreMx);
	Log("ES=%.4f\n", ES);
	FILE *f = CreateStdioFile(opt(scoremx));
	if (f != 0)
		fprintf(f, "// ES=%.4f\n", ES);
	FT.MxToSrc(f, FeatureName, ScoreMx);
	FT.MxToSrc(g_fLog, FeatureName, ScoreMx);
	CloseStdioFile(f);
	}