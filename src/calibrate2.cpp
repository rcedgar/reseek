#include "myutils.h"
#include "scop40bench.h"
#include "binner.h"

static const uint BinCount = 100;

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


	SCOP40Bench SB;
	SB.ReadChains(CalFN, "");

	DSSParams Params;
	Params.SetFromCmdLine(SB.GetDBSize());

	SB.Setup(Params);

	SB.m_QuerySelf = true;
	SB.m_ScoresAreEvalues = true;
	SB.Run();
	SB.SetTFs();
	const float MaxFPR = 0.1;
	SB.SetStats(MaxFPR, true);
	}
