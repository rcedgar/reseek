#include "myutils.h"
#include "calibratesearcher.h"
#include "chainreader2.h"
#include "binner.h"
#include "timing.h"
#include <set>

/***
* -calibrate3
* Calculate accumulated per-chain count bins for TS values.
* 32 bins, TS range -0.1 .. 0.3, written to -calib_output calibrate3.tsv
* Used by -calibrate4 to attempt similar linear fitting found to work
* well for all chains with supervised FPs. 
* Empirically, results are terrible.
***/

void CalibrateSearcher::WriteSlopeCalibOutput(FILE *f,
  uint BinCount, float TSlo, float TShi) const
	{
	if (f == 0)
		return;

// Mid-point TS values for bins
	vector<float> Mids;
	const vector<float> &TSs_0 = m_TestStatsVec[0];
	Binner<float> B(TSs_0, BinCount, TSlo, TShi);
	const vector<uint> &Bins = B.GetBins();
	//B.GetAccumBinsReverse(Bins);
	fprintf(f, "TS");
	for (uint Bin = 0; Bin < BinCount; ++Bin)
		fprintf(f, "\t%.3g", B.GetBinMid(Bin));
	fprintf(f, "\n");

	const uint NR = GetDBChainCount();
	for (uint i = 0; i < NR; ++i)
		{
		const char *Label = m_DBChains[i]->m_Label.c_str();
		const vector<float> &TSs_i = m_TestStatsVec[i];
		Binner<float> B(TSs_i, BinCount, TSlo, TShi);
		//B.GetAccumBinsReverse(Bins);
		const vector<uint> &Bins = B.GetBins();

		fprintf(f, "%s", Label);
		for (uint Bin = 0; Bin < BinCount; ++Bin)
			{
			uint n = Bins[Bin];
			fprintf(f, "\t%u", n);
			}
		fprintf(f, "\n");
		}
	}

/***
* Measure TS distribution.
* Unsupervised, includes TPs at high scores.
* Linear fit in calibrate4.
***/
void cmd_calibrate3()
	{
	if (optset_output)
		Die("Use -calib_output");

	optset_sensitive = true;
	opt_sensitive = true;

	optset_minchainlength = true;
	opt_minchainlength = 5;

	optset_evalue = true;
	opt_evalue = 10;

	const string &QFN = g_Arg1;
	const string &DBFN = g_Arg1;

	CalibrateSearcher DBS;
	DSSParams Params;
	Params.SetFromCmdLine(10000);
	DBS.m_Params = &Params;

	DBS.LoadDB(DBFN);
	DBS.Setup();

	ChainReader2 QCR;
	QCR.Open(QFN);
	DBS.RunQuery(QCR);

	DSSAligner::Stats();

	if (optset_calib_output)
		{
		FILE *fOut = CreateStdioFile(opt_calib_output);
		DBS.WriteSlopeCalibOutput(fOut, 33, 0, 0.3f);
		CloseStdioFile(fOut);
		}
	}
