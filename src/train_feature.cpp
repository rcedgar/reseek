#include "myutils.h"
#include "featuretrainer.h"

void cmd_train_feature()
	{
	asserta(optset_feature);
	asserta(optset_undef_binning);
	asserta(optset_output);
	asserta(!optset_alpha_sizes);

	optset_nofeatures = true;
	opt_nofeatures = true;

	const string FeatureName = string(opt(feature));
	FEATURE F = StrToFeature(FeatureName.c_str());

	const string &AlnsFN = g_Arg1;
	const string &ChainsFN = opt(db);

	const string &strUB = opt(undef_binning);
	UNDEF_BINNING UB = StrToUB(strUB);

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
	FT.m_UB = UB_IncludeUndefined;
	FT.SetInput(ChainsFN, AlnsFN);
	FT.SetFeature(F);

	uint BestDefaultLetter = UINT_MAX;
	if (UB == UB_UndefinedIsDefaultLetter)
		{
		FeatureTrainer FT2;
		FT2.SetOptionsFromCmdLine();
		FT2.m_UB = UB_IgnoreUndefined;
		FT2.SetInput(ChainsFN, AlnsFN);
		FT2.SetFeature(F);
		FT2.SetAlphaSize(AlphaSize, UB_IgnoreUndefined, UINT_MAX);
		FT2.Train();
		BestDefaultLetter = FT2.GetBestDefaultLetter(UINT_MAX);
		ProgressLog("BestDefaultLetter = %u\n",
					BestDefaultLetter);
		}

	FT.SetOptionsFromCmdLine();
	FT.SetAlphaSize(AlphaSize, UB, BestDefaultLetter);
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
	}
