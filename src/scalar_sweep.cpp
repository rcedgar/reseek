#include "myutils.h"
#include "scop40bench.h"
#include "features.h"

void cmd_scalar_sweep()
	{
	FILE *ftsv = CreateStdioFile(opt_output);
	asserta(optset_param);
	asserta(optset_minval);
	asserta(optset_maxval);
	asserta(optset_n);

	const string &CalFN = g_Arg1;
	SCOP40Bench SB;
	DSSParams Params;
	SB.Setup();
	SB.LoadDB(CalFN);

	uint Sens0 = 0;
	const uint StepCount = opt_n;
	float MinVal = (float) opt_minval;
	float MaxVal = (float) opt_maxval;
	for (uint Step = 0; Step <= StepCount; ++Step)
		{
		float Value = MinVal + Step*(MaxVal - MinVal)/StepCount;
		Progress("%s = %.4g\n", opt_param, Value);
		Params.SetParam(opt_param, Value, false);

		SB.RunSelf();
		uint Sens = SB.GetSens1stFP();

		char tofg = tof(opt_dpgaps);
		string Stem;
		GetStemName(CalFN, Stem);
		uint HitCount = SB.GetHitCount();
		ProgressLog("%s=%.3g\t%u\t%s\n", opt_param, Value, Sens, Stem.c_str());
		if (ftsv != 0)
			fprintf(ftsv, "%u\t%s\t%.4g\t%s\n", Sens, opt_param, Value, Stem.c_str());
		}
	CloseStdioFile(ftsv);
	}
