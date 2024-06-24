#include "myutils.h"
#include "calibratesearcher.h"

/***
* Calibrate distribution of FP errors on one_per_sf or one_per_fold.
* Gumbel fits well for the bulk of the distribution, but the fit is
* poor in both tails which causes E-value estimates to diverge in
* the most relevant ranges for practice.
***/
void cmd_calibrate()
	{
	const string &QCalFN = g_Arg1;
	const string &DBFN = opt_db;
	CalibrateSearcher DBS;
	DBS.LoadChains(QCalFN, DBFN);

	DSSParams Params;
	Params.SetFromCmdLine(DBS.GetDBSize());

	FILE *fOut = CreateStdioFile(opt_output);
	DBS.Setup(Params);
	DBS.Run();
	DSSAligner::Stats();
	DBS.ScanAll();
	DBS.SetAllBins();
	DBS.SetAllAccum();
	DBS.Setxys();
	DBS.FitGumbel();
	DBS.WriteBins(fOut);
	CloseStdioFile(fOut);
	}
