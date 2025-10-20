#include "myutils.h"
#include "featuretrainer2.h"
#include "sort.h"

void cmd_train_feature2()
	{
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

	const string &ChainFN = opt(db);			// "c:/src/reseek/test_data/scop40.bca";
	const string &TrainTPAlnFN = opt(traintps); // "../big_out/tp.a.mints05.maxts25.fa2";
	const string &EvalTPAlnFN = opt(evaltps);	// "../big_out/tp.a.evalrange.fa2";
	const string &EvalFPAlnFN = opt(evalfps);	// "../big_out/fp.a.evalrange.fa2";

///////////////////////////////////////////////////////////////////////////////////////
// Structure chains, must include all Train and Eval chains, may include others
///////////////////////////////////////////////////////////////////////////////////////
	vector<PDBChain *> Chains;
	map<string, uint> LabelToChainIdx;
	FeatureTrainer2::ReadChains(ChainFN, Chains, LabelToChainIdx);

///////////////////////////////////////////////////////////////////////////////////////
// Training alignments (TP only)
///////////////////////////////////////////////////////////////////////////////////////
	vector<bool> TrainsTPs_notused;
	vector<string> TrainRows;
	vector<string> TrainLabels;
	vector<uint> TrainChainIdxs;
	FeatureTrainer2::AppendAlns(TrainTPAlnFN, LabelToChainIdx, true,
	  TrainRows, TrainLabels, TrainChainIdxs, TrainsTPs_notused);

///////////////////////////////////////////////////////////////////////////////////////
// Eval alignments (TP and FP)
///////////////////////////////////////////////////////////////////////////////////////
	vector<bool> EvalTPs;
	vector<string> EvalRows;
	vector<string> EvalLabels;
	vector<uint> EvalRowChainIdxs;
	vector<uint> EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec;

	FeatureTrainer2::LoadEvalAlns(EvalTPAlnFN, EvalFPAlnFN, LabelToChainIdx,
		EvalRows, EvalLabels, EvalRowChainIdxs, EvalTPs,
		EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec);

	vector<vector<float> > ScoreMx;
	bool UndefsAllowed = true;
	float BestArea = FLT_MAX;
	vector<float> Areas;
	const string &BgMethod = opt(bgmethod);
	float BestOpenPenalty, BestExtPenalty, BestBias;
	if (opt(dss))
		{
		FeatureTrainer2::TrainDSSFeature(F, Chains, LabelToChainIdx,
			TrainRows, TrainLabels, TrainChainIdxs,
			EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
			EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
			BgMethod, ScoreMx,
			BestOpenPenalty, BestExtPenalty, BestBias, BestArea);
		return;
		}

	if (IsInt)
		{
		for (uint ReplaceUndefWithThisLetter = 0; ReplaceUndefWithThisLetter < AlphaSize;
			++ReplaceUndefWithThisLetter)
			{
			FeatureTrainer2::TrainIntFeature(F, Chains, LabelToChainIdx,
				TrainRows, TrainLabels, TrainChainIdxs,
				EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
				EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
				UndefsAllowed, ReplaceUndefWithThisLetter,
				BgMethod, ScoreMx, BestArea, 0);
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
		if (optset_rocsteps)
			{
			FILE *fSteps = CreateStdioFile(opt(rocsteps));
			uint BestLetter = Order[AlphaSize-1];
			FeatureTrainer2::TrainIntFeature(F, Chains, LabelToChainIdx,
				TrainRows, TrainLabels, TrainChainIdxs,
				EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
				EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
				UndefsAllowed, BestLetter,
				BgMethod, ScoreMx, BestArea, fSteps);
			CloseStdioFile(fSteps);
			}
		}
	else
		{
		FILE *fSteps = 0;
		if (optset_rocsteps)
			fSteps = CreateStdioFile(opt(rocsteps));

		float BestArea;
		FeatureTrainer2::TrainFloatFeature(
			F, AlphaSize, Chains, LabelToChainIdx,
			TrainRows, TrainLabels, TrainChainIdxs,
			EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
			EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
			ScoreMx, QS_UndefDistinct, BestArea, fSteps);
		CloseStdioFile(fSteps);
		}
	}
