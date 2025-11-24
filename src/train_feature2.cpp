#include "myutils.h"
#include "featuretrainer2.h"
#include "sort.h"

#define EVAL	0

void cmd_train_feature2()
	{
	const string FeatureName = g_Arg1;
	FEATURE F = StrToFeature(FeatureName.c_str());
	DSSParams::Init(DM_DefaultSensitive);

	const string &ChainFN = opt(db);			// "src/reseek/test_data/scop40.bca";
	const string &TrainTPAlnFN = opt(traintps); // "src/2025-10_reseek_tune/big_fa2/tp.mints05.maxts25.fa2";
	const string &EvalTPAlnFN = opt(evaltps);	// "src/2025-10_reseek_tune/big_fa2/tp.evalrange.fa2";
	const string &EvalFPAlnFN = opt(evalfps);	// "src/2025-10_reseek_tune/big_fa2/fp.evalrange.fa2";

	FILE *fOut = CreateStdioFile(opt(output));

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
	if (!opt(evalmu))
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

#if EVAL
	FeatureTrainer2::LoadEvalAlns(EvalTPAlnFN, EvalFPAlnFN, LabelToChainIdx,
		EvalRows, EvalLabels, EvalRowChainIdxs, EvalTPs,
		EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec);
#endif

	vector<vector<float> > ScoreMx;
	float BestArea = FLT_MAX;
	asserta(optset_evalmu || optset_background_style);
	BACKGROUND_STYLE BS = FeatureTrainer2::StrToBS(opt(background_style));
	FILE *fSteps = 0;
	if (opt(evalmu))
		{
		extern float ScoreMx_Mu[36][36];
		ScoreMx.resize(36);
		for (uint i = 0; i < 36; ++i)
			{
			ScoreMx[i].resize(36);
			for (uint j = 0; j < 36; ++j)
				ScoreMx[i][j] = ScoreMx_Mu[i][j];
			}
		uint ReplaceUndefWithThisLetter = 0; // see dss.cpp:697
		FeatureTrainer2::EvaluateMu(Chains, LabelToChainIdx,
			EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
			EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
			ScoreMx, ReplaceUndefWithThisLetter,
			BestArea, fSteps);
		}
	else if (opt(dss))
		{
		FeatureTrainer2::TrainDSSFeature(F, Chains, LabelToChainIdx,
			TrainRows, TrainLabels, TrainChainIdxs,
			EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
			EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
			ScoreMx, BS, 0);

		if (fOut != 0)
			{
			FeatureTrainer2::ScoreMxToTsv(fOut, ScoreMx);
			if (!FeatureIsInt(F))
				{
				vector<float> BinTs;
				DSSParams::GetBinTs(F, BinTs);
				FeatureTrainer2::BinTsToTsv(fOut, BinTs);
				}
			}
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
				BS, ScoreMx, 0, fSteps);
			FeatureTrainer2::ScoreMxToTsv(fOut, ScoreMx);
			FeatureTrainer2::ScoreMxToSrc(g_fLog, ScoreMx);
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

			vector<float> BinTs;
			FeatureTrainer2::TrainFloatFeature(
				F, AlphaSize, Chains, LabelToChainIdx,
				TrainRows, TrainLabels, TrainChainIdxs,
				EvalTPs, EvalRows, EvalLabels, EvalRowChainIdxs,
				EvalAlnColCountVec, EvalAlnOpenVec, EvalAlnExtVec,
				ScoreMx, BinTs, QS, UndefReplaceValue, 0, fSteps);
			FeatureTrainer2::ScoreMxToTsv(fOut, ScoreMx);
			FeatureTrainer2::BinTsToTsv(fOut, BinTs);
			FeatureTrainer2::BinTsToSrc(g_fLog, BinTs);
			FeatureTrainer2::ScoreMxToSrc(g_fLog, ScoreMx);
			}
		}
	Log("@FEV@ %s\n", FeatureTrainer2::m_FevStr.c_str());
	if (fOut != 0)
		{
		string CmdLine;
		GetCmdLine(CmdLine);
		fprintf(fOut, "cmd\t%s\n", CmdLine.c_str());
		fprintf(fOut, "fev\t%s\n", FeatureTrainer2::m_FevStr.c_str());
		fprintf(fOut, "git\t%s\n", g_GitVer);
		CloseStdioFile(fOut);
		}
	}
