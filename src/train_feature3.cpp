#include "myutils.h"
#include "letteralndb.h"
#include "pwalndb.h"
#include "bytevecdb.h"
#include "pdbchain.h"
#include "featuretrainer2.h"
#include "logmxptr.h"

void cmd_train_feature3()
	{
	const string FeatureName = g_Arg1;
	FEATURE F = StrToFeature(FeatureName.c_str());
	uint AS = DSSParams::GetAlphaSize(F);

	vector<PDBChain *> Chains;
	map<string, uint> LabelToChainIdx;
	FeatureTrainer2::ReadChains(opt(db), Chains, LabelToChainIdx);

	ByteVecDB BVDB;
	BVDB.Init(Chains, F);

	PWAlnDB PADB;
	PADB.Load(opt(input2), LabelToChainIdx, false);
	uint MissingCount = SIZE(PADB.m_MissingLabels);
	if (MissingCount > 0)
		ProgressLog("%u missing\n", MissingCount);

	LetterAlnDB LADB;
	LADB.Init(BVDB, PADB);
	LADB.SetCountsAndFreqs();

	LogVecPtr<float>("Freqs", "  %7.3g", LADB.m_FreqsPtr, AS);
	LogVecPtr<uint>("Counts", "  %7u", LADB.m_CountsPtr, AS);
	LogMxPtr<float>("JointFreqs", "  %7.3g", LADB.m_JointFreqsMxPtr, AS);
	LogMxPtr<uint>("JointCounts", "  %7u", LADB.m_JointCountsMxPtr, AS);

	uint PsuedoCount = 10;
	float *LOMx = LADB.GetLogOddsMxPtr(LADB.m_FreqsPtr, PsuedoCount);

	LogMxPtr<float>("LogOdds", "  %7.3g", LOMx, AS);
	}
