#include "myutils.h"
#include "features.h"
#include "pdbchain.h"
#include "dss.h"
#include "dssparams.h"
#include "binner.h"

const char *FeatureToStr(uint FeatureIndex)
	{
	return FeatureToStr(FEATURE(FeatureIndex));
	}

const char *FeatureToStr(FEATURE f)
	{
	switch (f)
		{
#define F(x)	case FEATURE_##x: return #x;
#include "featurelist.h"
		}
	Die("FeatureToStr(%u)", f);
	return "?";
	}

FEATURE StrToFeature(const char *s, bool ErrOk)
	{
#define F(x)	if (!stricmp(s, #x)) return FEATURE_##x;
#include "featurelist.h"
	if (!ErrOk)
		Die("StrToFeature(%s)", s);
	return FEATURE(UINT_MAX);
	}

uint StrToFeatureIndex(const char *s, bool ErrOk)
	{
#define F(x)	if (!stricmp(s, #x)) return uint(FEATURE_##x);
#include "featurelist.h"
	if (!ErrOk)
		Die("StrToFeatureIndex(%s)", s);
	return UINT_MAX;
	}

bool FeatureIsInt(FEATURE f)
	{
	switch (f)
		{
	case FEATURE_AA: return true;
	case FEATURE_B62: return true;

#define F(x)	case FEATURE_##x: return true;
#include "intfeatures.h"
#undef F

#define F(x)	case FEATURE_##x: return false;
#include "floatfeatures.h"
#undef F
		}
	asserta(false);
	return false;
	}

bool FeatureIsInt(uint FeatureIndex)
	{
	return FeatureIsInt(FEATURE(FeatureIndex));
	}

void cmd_feature_stats()
	{
	for (uint F = 0; F < FEATURE_COUNT; ++F)
		{
		uint AS = g_AlphaSizes2[F];
		float **Mx = g_ScoreMxs2[F];
		ProgressLog("[%2u]  %s",
		  F, FeatureToStr(F));
		if (Mx == 0)
			ProgressLog("  < missing scoremx");
		ProgressLog("\n");
		}
	}

void cmd_dump_float_feature()
	{
	const string FeatureName = g_Arg1;
	FEATURE F = StrToFeature(FeatureName.c_str());

	opt_force_undef = true;
	optset_force_undef = true;

	uint BinCount = 32;
	if (optset_bins)
		BinCount = opt(bins);

	vector<PDBChain *> Chains;
	ReadChains(opt(input), Chains);

	FILE *fOut = CreateStdioFile(opt(output));

	//DSSParams::Init(DM_UseCommandLineOption);
	DSS D;
	vector<float> Values;
	const uint ChainCount = SIZE(Chains);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *Chains[ChainIdx];
		D.Init(Chain);
		const uint L = Chain.GetSeqLength();
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			//float Value = D.GetFloatFeature(F, Pos);
			float ValueP = D.GetFloatFeature(FEATURE_PENDist, Pos);
			float ValueM = D.GetFloatFeature(FEATURE_MENDist, Pos);
			if (ValueP == FLT_MAX || ValueM == FLT_MAX)
				continue;
			float Value = ValueP - ValueM;
			if (Value == FLT_MAX)
				continue;
			Values.push_back(Value);
			if (fOut != 0)
				fprintf(fOut, "%.3g\n", Value);
			}
		}

	Binner<float> B(Values, BinCount);
	B.ToTsv(opt(output2));

	CloseStdioFile(fOut);
	}
