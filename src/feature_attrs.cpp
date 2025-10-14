#include "myutils.h"
#include "dssparams.h"
#include "dss.h"
#include "pdbchain.h"

void ReadChains(const string &FileName, vector<PDBChain *> &Chains);

static uint GetCounts(FEATURE F, vector<PDBChain *> &Chains,
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
			if (iForced == 0)
				Total += L;
			for (uint Pos = 0; Pos < L; ++Pos)
				{
				uint Letter = D.GetFeature(F, Pos);
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

	ProgressLog("%16.16s  %2.2s", "Feature", "AS");
	ProgressLog("  %5.5s     ", "Type");
	ProgressLog("  %18.18s", "Int forced");
	ProgressLog("  %18.18s", "Int not forced");
	ProgressLog("  %18.18s", "Float forced");
	ProgressLog("  %18.18s", "Float not forced");
	ProgressLog("\n");
	for (int iInt = 0; iInt < 2; ++iInt)
		{
		for (uint F = 0; F < FEATURE_COUNT; ++F)
			{
			if (iInt && !FeatureIsInt(F))
				continue;
			if (!iInt && FeatureIsInt(F))
				continue;
				
			uint AS = DSS::GetAlphaSize((FEATURE) F, true);
			if (AS == UINT_MAX)
				continue;
			uint UndefIntCountForced, UndefIntCountNotForced;
			uint UndefFloatCountForced, UndefFloatCountNotForced;
			uint Total = GetCounts((FEATURE) F, Chains,
				UndefIntCountForced, UndefIntCountNotForced,
				UndefFloatCountForced, UndefFloatCountNotForced);
			ProgressLog("%16.16s  %2u", FeatureToStr(F), AS);
			ProgressLog("  %5.5s", FeatureIsInt(F) ? "Int" : "Float");
			if (F == FEATURE_SS3 || F == FEATURE_NENSS3 || F == FEATURE_RENDist4)
				ProgressLog(" (Mu)");
			else
				ProgressLog("     ");

			ProgressLog("  %8u (%6.2f%%)", UndefIntCountForced, GetPct(UndefIntCountForced, Total));
			ProgressLog("  %8u (%6.2f%%)", UndefIntCountNotForced, GetPct(UndefIntCountNotForced, Total));
			if (!FeatureIsInt(F))
				{
				ProgressLog("  %8u (%6.2f%%)", UndefFloatCountForced, GetPct(UndefFloatCountForced, Total));
				ProgressLog("  %8u (%6.2f%%)", UndefFloatCountNotForced, GetPct(UndefFloatCountNotForced, Total));
				}
			ProgressLog("\n");
			}
		}
	}
