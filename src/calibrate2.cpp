#include "myutils.h"
#include "scop40bench.h"
#include "binner.h"

void cmd_calibrate2()
	{
	asserta(optset_benchlevel);

	string CalFN;
	if (g_Arg1 == ".")
#ifdef _MSC_VER
		CalFN = "c:/src/reseek_scop40/reseek_db/scop40_family.cal";
#else
		CalFN = "/c/src/reseek_scop40/reseek_db/scop40_family.cal";
#endif
	else
		CalFN = g_Arg1;


	const uint BinCount = 100;
	SCOP40Bench SB;
	SB.ReadChains(CalFN, "");

	DSSParams Params;
	Params.SetFromCmdLine(SB.GetDBSize());

	SB.Setup(Params);

	SB.m_QuerySelf = true;
	SB.m_ScoresAreEvalues = true;
	SB.Run();
	SB.SetTFs();

	const uint HitCount = SB.GetHitCount();
	asserta(SIZE(SB.m_TSs) == HitCount);
	asserta(SIZE(SB.m_TFs) == HitCount);
	vector<float> MinusLogTSs;
	vector<float> TSs;
	for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
		{
		if (SB.m_TFs[HitIdx])
			continue;
		float TS = SB.m_TSs[HitIdx];
		asserta(TS > 0);
		TSs.push_back(TS);
		float MinusLogTS = -logf(TS);
		MinusLogTSs.push_back(MinusLogTS);
		}

	Binner<float> Blin(TSs, BinCount, 0);
	Binner<float> Blog(MinusLogTSs, BinCount, 0);
	Blin.ToTsv(g_fLog);
	}
