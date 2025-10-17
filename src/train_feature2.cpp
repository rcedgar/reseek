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
	const string &TrainTPAlnFN = opt(input); // "../big_out/tp.a.mints05.maxts25.fa2";
	const string &TrainFPAlnFN = opt(input2); // "../big_out/fp.a.mints05.maxts25.fa2";
	const string &EvalTPAlnFN = TrainTPAlnFN;
	const string &EvalFPAlnFN = TrainFPAlnFN;
	vector<vector<float> > ScoreMx;
	FeatureTrainer2::SetIntFeature(F);
	bool UndefsAllowed = true;
	vector<float> Areas;
	const string &BgMethod = opt(bgmethod);
	for (uint ReplaceUndefWithThisLetter = 0; ReplaceUndefWithThisLetter < 16;
		++ReplaceUndefWithThisLetter)
		{
		float BestArea;
		FeatureTrainer2::TrainIntFeature(F, ChainFN,
			TrainTPAlnFN, EvalTPAlnFN, EvalFPAlnFN,
			UndefsAllowed, ReplaceUndefWithThisLetter,
			BgMethod, ScoreMx, BestArea);
		Areas.push_back(BestArea);
		}

	Log("\n");
	vector<uint> Order(AlphaSize);
	QuickSortOrder(Areas.data(), AlphaSize, Order.data());
	for (uint k = 0; k < AlphaSize; ++k)
		{
		uint Letter = Order[k];
		Log("[%2u]  %6.4f\n", Letter, Areas[Letter]);
		}
	}
