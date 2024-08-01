#include "myutils.h"
#include "dssparams.h"

double GetPNN(uint Idx, double TS);

void cmd_evalue_report()
	{
	DSSParams Params;
	Params.SetFromCmdLine(10000);

	const float LoTS = 0;
	const float HiTS = 1.55f;
	const float dTS = 0.1f;
	const float m = 20.5f;
	const float b = 2.9f;

	Log("\n");
	float TS = LoTS;
	for (;;)
		{
		float PredMinusLogP_error = m*TS + b;
		float P_error = expf(-PredMinusLogP_error);
		float E = Params.GetEvalue(TS);


		Log("TS = %4.2f", TS);
		Log("  P_error %8.3g", P_error);
		Log("  E %8.3g", E);

		for (uint Idx = 0; Idx < 4; ++Idx)
			{
			double PNN = GetPNN(Idx, TS);
			Log("  P_nn %8.3g", 1-PNN);
			}

		for (uint Idx = 0; Idx < 4; ++Idx)
			{
			double PNN = GetPNN(Idx, TS);
			Log("  Enn %8.3g", (1-PNN)*10000);
			}
		Log("\n");

		TS += dTS;
		if (TS > HiTS)
			break;
		}
	}
