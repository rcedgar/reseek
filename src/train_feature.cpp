#include "myutils.h"
#include "featuretrainer.h"

void cmd_train_feature()
	{
	asserta(optset_feature);
	asserta(optset_undef_binning);
	asserta(optset_output);

	opt_force_undef_letters = true;
	optset_force_undef_letters = true;

	const string FeatureName = string(opt(feature));
	FEATURE F = StrToFeature(FeatureName.c_str());

	const string &AlnsFN = g_Arg1;
	const string &ChainsFN = opt(db);

	const string &strUB = opt(undef_binning);
	UNDEF_BINNING UB = StrToUB(strUB);

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
	FT.SetInput(ChainsFN, AlnsFN);
	FT.SetFeature(F);
	if (!FeatureIsInt(F))
		{
		uint UndefCount = 0;
		FT.SetFloatValues(FLT_MAX, UndefCount);
		uint BinCount = FT.GetBinThresholdCount(AlphaSize, UB_NeverUndefined);
		FT.CalcBinTs(BinCount);
		}

#if	0
	uint BestDefaultLetter = UINT_MAX;
	if (FT.MustOptimizeDefaultLetter(UB))
		{
		uint UndefCount = 0;
		FT.SetFloatValues(FLT_MAX, UndefCount);
		uint BinCount = FT.GetBinThresholdCount(AlphaSize, UB_NeverUndefined);
		FT.CalcBinTs(BinCount);
		BestDefaultLetter = FT.GetBestDefaultLetter(UINT_MAX);
		}

	FT.Train();
	FT.WriteSummary(g_fLog);
	FT.WriteSummary(stderr);
	FT.ToSrc(g_fLog);

	FT.ToTsv(opt(output));

	vector<float> ExpectedScores;
	FT.GetExpectedScores(ExpectedScores);
	Log("\n");
	Log("%s(%u) expected scores\n", FeatureName.c_str(), AlphaSize);
	asserta(SIZE(ExpectedScores) == AlphaSize);
	for (uint i = 0; i < AlphaSize; ++i)
		Log("%2u = %.3g\n", i, ExpectedScores[i]);
#endif
	}
