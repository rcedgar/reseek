#include "myutils.h"
#include "featuretrainer.h"

void cmd_train_feature()
	{
	asserta(optset_feature);
//	asserta(optset_output);

	const string FeatureName = string(opt(feature));
	FEATURE F = StrToFeature(FeatureName.c_str());

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
	FT.SetInput(ChainsFN, AlnsFN);
	FT.SetFeature(F, AlphaSize);
	if (!FeatureIsInt(F))
		FT.TrainQuantization();
	FT.TrainLogOdds(false);
	}
