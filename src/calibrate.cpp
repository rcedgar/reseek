#include "myutils.h"
#include "chainreader2.h"
#include "calibratesearcher.h"
#include "statsig.h"

/***
* Calibrate distribution of FP errors on one_per_sf or one_per_fold.
* Gumbel fits well for the bulk of the distribution, but the fit is
* poor in both tails which causes E-value estimates to diverge in
* the most relevant ranges for practice.
***/
void cmd_calibrate()
	{
	optset_sensitive = true;
	opt(sensitive) = true;

	optset_minchainlength = true;
	opt(minchainlength) = 1;

	const string &QFN = g_Arg1;
	const string &DBFN = g_Arg1;

	CalibrateSearcher DBS;
	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast);
	DBS.m_Params = &Params;

	DBS.LoadDB(DBFN);
	DBS.Setup();

	ChainReader2 QCR;
	QCR.Open(QFN);
	DBS.RunQuery(QCR);

	FILE *fOut = CreateStdioFile(opt(output));
	DSSAligner::Stats();
	DBS.ScanAll();
	DBS.SetAllBins();
	DBS.SetAllAccum();
	DBS.Setxys();
	DBS.FitGumbel();
	DBS.WriteBins(fOut);
	CloseStdioFile(fOut);
	}
