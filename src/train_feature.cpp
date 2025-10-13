#include "myutils.h"
#include "featuretrainer.h"

void cmd_train_feature()
	{
	asserta(optset_feature);
	DSSParams::Init(DM_DefaultFast);

	opt_force_undef = true;
	optset_force_undef = true;

	const string FeatureName = string(opt(feature));
	FEATURE F = StrToFeature(FeatureName.c_str());
	asserta(optset_unaligned_background);
	asserta(optset_undef_overlap);

	bool UseUnalignedBackground = false;
	if (string(opt(unaligned_background)) == "yes")
		UseUnalignedBackground = true;
	else if (string(opt(unaligned_background)) == "no")
		UseUnalignedBackground = false;
	else
		Die("Invalid -unaligned_background %s", opt(unaligned_background));

	bool UndefOverlap = false;
	if (string(opt(undef_overlap)) == "yes")
		UndefOverlap = true;
	else if (string(opt(undef_overlap)) == "no")
		UndefOverlap = false;
	else
		Die("Invalid -undef_overlap %s", opt(undef_overlap));

	DSSParams::Init(DM_DefaultSensitive);

	uint AlphaSize = UINT_MAX;
	bool IsInt = FeatureIsInt(F);
	if (IsInt)
		AlphaSize = DSS::GetAlphaSize(F);
	else
		{
		asserta(optset_alpha_size);
		AlphaSize = opt(alpha_size);
		}

	FeatureTrainer FT;
	FT.m_UseUnalignedBackground = UseUnalignedBackground;
	FT.ReadChains(g_Arg1);
	FT.ReadAlns(opt(input), true);
	if (optset_input2)
		FT.ReadAlns(opt(input2), false);
	ProgressLog("%u TPs, %u FPs\n", FT.m_TPCount, FT.m_FPCount);
	FT.SetFeature(F, AlphaSize);
	if (FeatureIsInt(F))
		{
		if (UndefOverlap)
			FT.TrainInt_UndefOverlap();
		else
			FT.TrainInt_UndefDistinct();
		}
	else
		{
		if (UndefOverlap)
			FT.TrainFloat_UndefOverlap();
		else
			FT.TrainFloat_UndefDistinct();
		}
	FT.SetAlnSubstScores();

	float MinGapOpen = -0.2;
	float MinGapExt = -0.2;
	float dGapOpen = 0.2;
	float dGapExt = 0.05;
	float BestArea = 0;
	float BestGapOpen = 0;
	float BestGapExt = 0;
	for (uint i = 0; i < 10; ++i)
		for (uint j = 0; j < 10; ++j)
			{
			float GapOpen = MinGapOpen + i*dGapOpen;
			float GapExt = MinGapExt + j*dGapExt;
			FT.SetAlnScoresAndArea(GapOpen, GapExt);
			if (FT.m_Area > BestArea)
				{
				BestArea = FT.m_Area;
				BestGapOpen = GapOpen;
				BestGapExt = GapExt;
				}
			}
	ProgressLog("Best area %.4f, open %.4g, ext %.4g\n",
		BestArea, BestGapOpen, BestGapExt);

	FT.OptimizeGapPenalties();

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