#include "myutils.h"
#include "letteralndb.h"
#include "pwalndb.h"
#include "bytevecdb.h"
#include "pdbchain.h"
#include "featuretrainer2.h"
#include "logmxptr.h"

static void TrainTPFP()
	{
	asserta(optset_fasta2_tp);
	asserta(optset_fasta2_fp);
	const string FeatureName = g_Arg1;
	FEATURE F = StrToFeature(FeatureName.c_str());
	uint AS = DSSParams::GetAlphaSize(F);

	vector<PDBChain *> Chains;
	map<string, uint> LabelToChainIdx;
	FeatureTrainer2::ReadChains(opt(db), Chains, LabelToChainIdx);

	ByteVecDB BVDB;
	BVDB.Init(Chains, F);

	PWAlnDB PA_TP, PA_FP;
	PA_TP.Load(opt(fasta2_tp), LabelToChainIdx, false);
	PA_FP.Load(opt(fasta2_fp), LabelToChainIdx, false);
	uint MissingCount_TP = SIZE(PA_TP.m_MissingLabels);
	uint MissingCount_FP = SIZE(PA_FP.m_MissingLabels);
	if (MissingCount_TP > 0 || MissingCount_FP > 0)
		ProgressLog("%u TPs %u FPs missing\n",
			MissingCount_TP, MissingCount_FP);

	uint PsuedoCount = 10;
	LetterAlnDB LA_TP, LA_FP;
	LA_TP.Init(BVDB, PA_TP);
	LA_FP.Init(BVDB, PA_FP);
	LA_TP.SetCountsAndFreqs(PsuedoCount);
	LA_FP.SetCountsAndFreqs(PsuedoCount);
	ProgressLog("TP letter pairs %s bytes\n", MemBytesToStr(LA_TP.m_LetterPairCount*2));
	ProgressLog("FP letter pairs %s bytes\n", MemBytesToStr(LA_FP.m_LetterPairCount*2));

	LogVecPtr<float>("Freqs_TP", "  %7.3g", LA_TP.m_FreqsPtr, AS);
	LogVecPtr<uint>("Counts_TP", "  %7u", LA_TP.m_CountsPtr, AS);
	LogMxPtr<float>("JointFreqs_TP", "  %7.3g", LA_TP.m_JointFreqsMxPtr, AS);
	LogMxPtr<uint>("JointCounts_TP", "  %7u", LA_TP.m_JointCountsMxPtr, AS);

	LogVecPtr<float>("Freqs_FP", "  %7.3g", LA_FP.m_FreqsPtr, AS);
	LogVecPtr<uint>("Counts_FP", "  %7u", LA_FP.m_CountsPtr, AS);
	LogMxPtr<float>("JointFreqs_FP", "  %7.3g", LA_FP.m_JointFreqsMxPtr, AS);
	LogMxPtr<uint>("JointCounts_FP", "  %7u", LA_FP.m_JointCountsMxPtr, AS);

	float *LO_TP = LA_TP.GetLogOddsMxPtr(LA_TP.m_FreqsPtr);
	LogMxPtr<float>("LogOdds_TP", "  %7.3g", LO_TP, AS);

	float ES_TP = LetterAlnDB::GetES2(AS, LO_TP, LA_TP.m_JointFreqsMxPtr);
	float ES_FP = LetterAlnDB::GetES2(AS, LO_TP, LA_FP.m_JointFreqsMxPtr);
	Log("ES_TP=%.3g ES_FP=%.3g\n", ES_TP, ES_FP);
	}

void cmd_train_feature3()
	{
	asserta(optset_fasta2_tp);
	if (optset_fasta2_fp)
		{
		TrainTPFP();
		return;
		}

	const string FeatureName = g_Arg1;
	FEATURE F = StrToFeature(FeatureName.c_str());
	uint AS = DSSParams::GetAlphaSize(F);

	vector<PDBChain *> Chains;
	map<string, uint> LabelToChainIdx;
	FeatureTrainer2::ReadChains(opt(db), Chains, LabelToChainIdx);

	ByteVecDB BVDB;
	BVDB.Init(Chains, F);

	PWAlnDB PADB;
	PADB.Load(opt(fasta2_tp), LabelToChainIdx, false);
	uint MissingCount = SIZE(PADB.m_MissingLabels);
	if (MissingCount > 0)
		ProgressLog("%u missing\n", MissingCount);

	uint PsuedoCount = 10;
	LetterAlnDB LADB;
	LADB.Init(BVDB, PADB);
	LADB.SetCountsAndFreqs(PsuedoCount);
	ProgressLog("Letter pairs %s bytes\n", MemBytesToStr(LADB.m_LetterPairCount*2));

	LogVecPtr<float>("Freqs", "  %7.3g", LADB.m_FreqsPtr, AS);
	LogVecPtr<uint>("Counts", "  %7u", LADB.m_CountsPtr, AS);
	LogMxPtr<float>("JointFreqs", "  %7.3g", LADB.m_JointFreqsMxPtr, AS);
	LogMxPtr<uint>("JointCounts", "  %7u", LADB.m_JointCountsMxPtr, AS);

	float *LOMx = LADB.GetLogOddsMxPtr(LADB.m_FreqsPtr);

	LogMxPtr<float>("LogOdds", "  %7.3g", LOMx, AS);
	ProgressLog("ES %8.3f  %s\n", LADB.GetES(LADB.m_FreqsPtr), FeatureToStr(F));
	}
