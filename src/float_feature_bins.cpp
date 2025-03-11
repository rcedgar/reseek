#include "myutils.h"
#include "pdbchain.h"
#include "seqdb.h"
#include "alpha.h"
#include "dss.h"
#include "logodds.h"
#include "scop40bench.h"

static uint g_FeatureIndex = UINT_MAX;

static void GetBins(const vector<float> &Values, uint AlphaSize,
					vector<float> &BinTs)
	{
	const uint K = SIZE(Values);
	asserta(K > 0);
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		{
		uint k = ((i+1)*K)/AlphaSize;
		float t = Values[k];
		BinTs.push_back(t);
		}
	asserta(SIZE(BinTs) == AlphaSize - 1);
	}

static void ReportBins(const vector<float> &Values, uint AlphaSize)
	{
	const char *FeatureName = FeatureToStr(g_FeatureIndex);
	vector<float> BinTs;
	GetBins(Values, AlphaSize, BinTs);

	ProgressLog("%s: AlphaSize %u\n", FeatureName, AlphaSize);

	Log("\n// %s [%2u]\n", FeatureName, AlphaSize);
	Log("ALPHA_SIZE(%s, %u);\n", FeatureName, AlphaSize);

	Log("BIN_T_BEGIN(%s);\n", FeatureName);
	for (uint i = 0; i + 1 < AlphaSize; ++i)
		Log("BIN_T(%s, %u, %.4g);\n",
		  FeatureName, i, BinTs[i]);
	Log("BIN_T_END(%s);\n", FeatureName);
	}

void GetFloatFeatureValues(const vector<PDBChain *> &Chains, 
						   FEATURE F, vector<float> &Values)
	{
	Values.clear();

	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast, SCOP40_DBSIZE);

	DSS D;
	D.SetParams(Params);

	const uint ChainCount = SIZE(Chains);
	for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, ChainCount, "Processing");
		const PDBChain &Chain = *Chains[ChainIndex];
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			float Value = D.GetFloatFeature(F, Pos);
			Values.push_back(Value);
			}
		}
	const uint K = SIZE(Values);
	asserta(K > 0);
	sort(Values.begin(), Values.end());
	}

void ConstructBinsFromValues(const vector<float> &Values, uint AlphaSize,
							 vector<float> &BinTs)
	{
	GetBins(Values, AlphaSize, BinTs);
	}

void cmd_float_feature_bins()
	{
	asserta(optset_feature);
	g_FeatureIndex = StrToFeatureIndex(opt(feature));
	const char *FeatureName = FeatureToStr(g_FeatureIndex);
	FEATURE F = FEATURE(g_FeatureIndex);

	vector<PDBChain *> Chains;
	ReadChains(opt(train_cal), Chains);
	const uint ChainCount = SIZE(Chains);

	vector<float> Values;
	GetFloatFeatureValues(Chains, F, Values);

	if (optset_alpha_size)
		ReportBins(Values, opt(alpha_size));
	else
		{
		ReportBins(Values, 3);
		ReportBins(Values, 4);
		ReportBins(Values, 6);
		ReportBins(Values, 8);
		ReportBins(Values, 10);
		ReportBins(Values, 12);
		ReportBins(Values, 16);
		ReportBins(Values, 24);
		ReportBins(Values, 32);
		}
	}
