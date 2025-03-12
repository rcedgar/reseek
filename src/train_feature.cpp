#include "myutils.h"
#include "featuretrainer.h"

void cmd_train_feature()
	{
	asserta(optset_feature);
	const string FeatureName = string(opt(feature));
	FEATURE F = StrToFeature(FeatureName.c_str());
	FILE *fOut = CreateStdioFile(opt(output));

	const string &AlnsFN = g_Arg1;
	const string &ChainsFN = opt(db);

	vector<uint> AlphaSizes;
	bool IsInt = FeatureIsInt(F);
	if (IsInt)
		{
		uint AS = GetAlphaSize(F);
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
	FT.Init(ChainsFN, AlnsFN);
	FT.SetFeature(F);

	const uint ASCount = SIZE(AlphaSizes);
	for(uint i = 0; i < ASCount; ++i)
		{
		uint AS = AlphaSizes[i];
		FT.SetAlphaSize(AS);
		FT.Train();
		FT.WriteSummary(g_fLog);
		FT.WriteSummary(stderr);
		if (i == 0)
			FT.WriteTsvHdr(fOut, ASCount);
		FT.ToTsv(fOut);
		}
	CloseStdioFile(fOut);
	}
