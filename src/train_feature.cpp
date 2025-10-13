#include "myutils.h"
#include "featuretrainer.h"

void cmd_train_feature()
	{
	asserta(optset_feature);
	DSSParams::Init(DM_DefaultFast);

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

	const string &AlnsFN = g_Arg1;
	const string &ChainsFN = opt(db);

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
	FT.SetInput(ChainsFN, AlnsFN);
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
	FT.ToTsv(opt(output));
	vector<vector<float> > ScoreMx;
	float ES = FT.GetLogOddsMx(ScoreMx);
	FILE *f = CreateStdioFile(opt(scoremx));
	if (f != 0)
		fprintf(f, "// ES=%.4f\n", ES);
	FT.MxToSrc(f, FeatureName, ScoreMx);
	CloseStdioFile(f);
	}
