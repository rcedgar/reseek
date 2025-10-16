#include "myutils.h"
#include "featuretrainer2.h"
#include "sort.h"

void cmd_train_feature2()
	{
	DSSParams::Init(DM_DefaultFast);

	const string FeatureName = g_Arg1;
	FEATURE F = StrToFeature(FeatureName.c_str());

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

	const string &ChainFN = "c:/src/reseek/test_data/scop40.bca";
	const string &TrainTPAlnFN = "../big_out/tp.a.mints05.maxts25.fa2";
	const string &TrainFPAlnFN = "../big_out/tp.b.mints05.maxts25.fa2";
	const string &EvalTPAlnFN = TrainTPAlnFN;
	const string &EvalFPAlnFN = TrainFPAlnFN;
	vector<vector<float> > ScoreMx;
	FeatureTrainer2::SetIntFeature(F);
	FeatureTrainer2::TrainIntFeatureNoUndefs(F, ChainFN,
		TrainTPAlnFN, TrainFPAlnFN, EvalTPAlnFN, EvalFPAlnFN,
		ScoreMx);
	}
