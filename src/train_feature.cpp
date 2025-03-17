#include "myutils.h"
#include "featuretrainer.h"

void cmd_train_feature()
	{
	asserta(optset_feature);
	asserta(optset_wildcard);
	asserta(!optset_output);

	const string FeatureName = string(opt(feature));
	FEATURE F = StrToFeature(FeatureName.c_str());

	const string &AlnsFN = g_Arg1;
	const string &ChainsFN = opt(db);

	const string strWildcard = opt(wildcard);
	bool Wildcard = false;
	if (strWildcard == "yes")
		Wildcard = true;
	else if (strWildcard == "no")
		Wildcard = false;
	else
		Die("Bad -wildcard");

	vector<uint> AlphaSizes;
	bool IsInt = FeatureIsInt(F);
	if (IsInt)
		{
		uint AS = DSS::GetAlphaSize(F);
		asserta(AS != 0 && AS != UINT_MAX);
		AlphaSizes.push_back(AS);
		}
	else
		{
		if (optset_alpha_size)
			{
			uint AS = opt(alpha_size);
			asserta(AS != 0 && AS != UINT_MAX);
			AlphaSizes.push_back(AS);
			}
		else if (optset_alpha_sizes)
			{
			vector<string> Fields;
			Split(opt(alpha_sizes), Fields, ',');
			for (uint i = 0; i < SIZE(Fields); ++i)
				{
				uint AS = StrToUint(Fields[i]);
				AlphaSizes.push_back(AS);
				}
			}
		else
			Die("Must set -alpha_size[s]");
		}
	asserta(!AlphaSizes.empty());

	FeatureTrainer FT;
	FT.SetInput(ChainsFN, AlnsFN);
	FT.SetFeature(F);

	const uint ASCount = SIZE(AlphaSizes);
	for(uint i = 0; i < ASCount; ++i)
		{
		uint AS = AlphaSizes[i];
		FT.SetAlphaSize(AS);
		FT.Train(Wildcard);
		FT.WriteSummary(g_fLog);
		FT.WriteSummary(stderr);

		string OutputFN;
		Ps(OutputFN, "%s.%u", FeatureName.c_str(), AS);
		FT.ToTsv(OutputFN);

		vector<float> ExpectedScores;
		FT.GetExpectedScores(ExpectedScores);
		Log("\n");
		Log("%s(%u) expected scores\n", FeatureName.c_str(), AS);
		asserta(SIZE(ExpectedScores) == AS);
		for (uint i = 0; i < AS; ++i)
			Log("%2u = %.3g\n", i, ExpectedScores[i]);
		}
	}
