#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "pdbchain.h"

void ReadChains(const string &FileName, vector<PDBChain *> &Chains);

static uint GetCounts(FEATURE F, vector<PDBChain *> &Chains,
	vector<uint> &Counts,
	uint &UndefIntCountForced, uint &UndefIntCountNotForced,
	uint &UndefFloatCountForced, uint &UndefFloatCountNotForced)
	{
	uint Total = 0;
	UndefIntCountForced = 0;
	UndefIntCountNotForced = 0;
	UndefFloatCountForced = 0;
	UndefFloatCountNotForced = 0;

	uint AS = DSS::GetAlphaSize(F);
	bool IsInt = FeatureIsInt(F);

	DSS D;
	const uint ChainCount = SIZE(Chains);
	Counts.clear();
	Counts.resize(AS);
	for (int iForced = 0; iForced < 2; ++iForced)
		{
		if (iForced)
			{
			opt_force_undef = true;
			optset_force_undef = true;
			}
		else
			{
			opt_force_undef = false;
			optset_force_undef = false;
			}
		uint UndefCount = 0;
		for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
			{
			const PDBChain &Chain = *Chains[ChainIdx];
			D.Init(Chain);
			const uint L = Chain.GetSeqLength();
			Total += L;
			for (uint Pos = 0; Pos < L; ++Pos)
				{
				uint Letter = D.GetFeature(F, Pos);
				if (Letter < AS)
					Counts[Letter] += 1;
				if (Letter == UINT_MAX)
					{
					if (iForced == 1)
						++UndefIntCountForced;
					else
						++UndefIntCountNotForced;
					}
				if (!IsInt)
					{
					float Value = D.GetFloatFeature(F, Pos);
					if (Value == FLT_MAX)
						{
						if (iForced == 1)
							++UndefFloatCountForced;
						else
							++UndefFloatCountNotForced;
						}
					}
				}
			}
		}
	return Total;
	}

void cmd_feature_attrs()
	{
	DSSParams::Init(DM_DefaultFast);

	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);

	for (uint F = 0; F < FEATURE_COUNT; ++F)
		{
		uint AS = DSS::GetAlphaSize((FEATURE) F, true);
		if (AS == UINT_MAX)
			continue;
		vector<uint> Counts;
		uint UndefIntCountForced, UndefIntCountNotForced;
		uint UndefFloatCountForced, UndefFloatCountNotForced;
		uint Total = GetCounts((FEATURE) F, Chains, Counts,
			UndefIntCountForced, UndefIntCountNotForced,
			UndefFloatCountForced, UndefFloatCountNotForced);
		asserta(SIZE(Counts) == AS);
		ProgressLog("%16.16s  %2u", FeatureToStr(F), AS);
		ProgressLog("  %10u", Total);
		ProgressLog("  %10u", UndefIntCountForced);
		ProgressLog("  %10u", UndefIntCountNotForced);
		ProgressLog("  %10u", UndefFloatCountForced);
		ProgressLog("  %10u", UndefFloatCountNotForced);
		ProgressLog("\n");
		}
	}
