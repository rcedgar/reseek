#include "myutils.h"
#include "featuretrainer2.h"
#include "sort.h"

void cmd_train_feature2()
	{
	const string FeatureName = g_Arg1;
	FEATURE F = StrToFeature(FeatureName.c_str());
	DSSParams::Init(DM_DefaultSensitive);

	const string &ChainFN = opt(db);			// "c:/src/reseek/test_data/scop40.bca";
	const string &TrainTPAlnFN = opt(traintps); // "../big_out/tp.a.mints05.maxts25.fa2";
	const string &EvalTPAlnFN = opt(evaltps);	// "../big_out/tp.a.evalrange.fa2";
	const string &EvalFPAlnFN = opt(evalfps);	// "../big_out/fp.a.evalrange.fa2";

	FeatureTrainer2::m_FevStr.clear();

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
	FeatureTrainer2::AppendAlns("traintps", TrainTPAlnFN, LabelToChainIdx, true,
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
	float BestArea = FLT_MAX;
	asserta(optset_background_style);
	BACKGROUND_STYLE BS = FeatureTrainer2::StrToBS(opt(background_style));
	FILE *fSteps = 0;
	if (opt(dss))
		{
		FeatureTrainer2::TrainDSSFeature(F, Chains, LabelToChainIdx,
			TrainRows, TrainLabels, TrainChainIdxs,
			EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
			EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
			ScoreMx, BS, BestArea);
		}
	else
		{
		uint AlphaSize = UINT_MAX;
		bool IsInt = FeatureIsInt(F);
		if (IsInt)
			AlphaSize = DSSParams::GetAlphaSize(F);
		else
			{
			asserta(optset_alpha_size);
			AlphaSize = opt(alpha_size);
			}

		if (IsInt)
			{
			asserta(optset_undef_letter);
			uint ReplaceUndefWithThisLetter = opt(undef_letter);
			bool UndefsAllowed = true;
			FeatureTrainer2::TrainIntFeature(F, Chains, LabelToChainIdx,
				TrainRows, TrainLabels, TrainChainIdxs,
				EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
				EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
				UndefsAllowed, ReplaceUndefWithThisLetter,
				BS, ScoreMx, BestArea, fSteps);
			}
		else
			{
			asserta(BS == BS_Float);
				asserta(optset_quantize_style);
			float UndefReplaceValue = FLT_MAX;
			QUANTIZE_STYLE QS = FeatureTrainer2::StrToQS(opt(quantize_style));
			if (QS == QS_UndefReplaceUser)
				{
				asserta(optset_undef_value);
				UndefReplaceValue = (float) opt(undef_value);
				}
			FeatureTrainer2::TrainFloatFeature(
				F, AlphaSize, Chains, LabelToChainIdx,
				TrainRows, TrainLabels, TrainChainIdxs,
				EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
				EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
				ScoreMx, QS, UndefReplaceValue, BestArea, fSteps);
			}
		}
	Log("@FEV@ %s\n", FeatureTrainer2::m_FevStr.c_str());
	}
