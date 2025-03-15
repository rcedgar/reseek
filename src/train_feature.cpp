#include "myutils.h"
#include "featuretrainer.h"

void cmd_train_feature()
	{
	asserta(optset_feature);

	const string FeatureName = string(opt(feature));
	FEATURE F = StrToFeature(FeatureName.c_str());

	const string &AlnsFN = g_Arg1;
	const string &ChainsFN = opt(db);

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
		FT.Train();
		FT.WriteSummary(g_fLog);
		FT.WriteSummary(stderr);

		string OutputFN;
		Ps(OutputFN, "%s.%u", FeatureName.c_str(), AS);
		if (optset_output)
			OutputFN = string(opt(output)) + "." + OutputFN;
		FT.ToTsv(OutputFN);
		}
	}
