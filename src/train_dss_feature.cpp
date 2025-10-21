#include "myutils.h"
#include "featuretrainer2.h"
#include "sort.h"

void cmd_train_dss_feature()
	{
	Die("Obsolete use train_feature2 -dss");
#if 0
	const string FeatureName = g_Arg1;
	FEATURE F = StrToFeature(FeatureName.c_str());
	DSSParams::Init(DM_DefaultSensitive);

	const string &ChainFN = opt(db); // "c:/src/reseek/test_data/scop40.bca";
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

	uint AS = DSS::GetAlphaSize(F, true);
	if (AS == UINT_MAX)
		Die("Not DSS feature %s", FeatureName.c_str());

	if (FeatureIsInt(F))
		FeatureTrainer2::SetIntFeature(F);
	else
		{
		uint AlphaSize = DSS::GetAlphaSize(F);
		FeatureTrainer2::SetFloatFeature(F, AlphaSize);
		}

	float BestOpenPenalty;
	float BestExtPenalty;
	float BestBias;
	float BestArea;
	vector<vector<float> > ScoreMx;

	extern float **g_ScoreMxs2[FEATURE_COUNT];
	const float * const *Mx = g_ScoreMxs2[F];
	if (F == FEATURE_Mu || Mx != 0)
		{
		FeatureTrainer2::m_BS = BS_DSSScoreMx;
		vector<vector<uint> > ChainIntSeqsNoUndefs;
		FeatureTrainer2::GetDSSScoreMx(F, ScoreMx);
		FeatureTrainer2::GetChainIntSeqs_DSS(Chains, ChainIntSeqsNoUndefs);
		FeatureTrainer2::LogChainIntSeqsStats(ChainIntSeqsNoUndefs);
		FeatureTrainer2::EvalLogOddsMx(ChainIntSeqsNoUndefs, EvalRows, EvalRowChainIdxs,
			EvalTPs, EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
			ScoreMx, BestOpenPenalty, BestExtPenalty, BestBias, BestArea);
		}

	for (uint BgMethodIdx = 0; BgMethodIdx < 2; ++BgMethodIdx)
		{
		BACKGROUND_STYLE BS = BS_Invalid;
		if (BgMethodIdx == 0)
			BS = BS_AlignedLetters;
		else if (BgMethodIdx == 1)
			BS = BS_UniqueChains;
		else
			asserta(false);
		FeatureTrainer2::TrainDSSFeature(F, Chains, LabelToChainIdx,
			TrainRows, TrainLabels, TrainChainIdxs,
			EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
			EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
			ScoreMx, BS, BestArea);
		}
#endif
	}
