#include "myutils.h"
#include "scop40bench.h"
#include "binner.h"
#include "timing.h"
#include <set>

/***
* -calibrate5
* Per-chain fit log(TS) to Gumbel
* y is PDF after trimming outliers (putative TPs)
***/

#define LOG_FIT	0
static const uint NOUTLIERS = 3;
static const uint NBINS = 32;

double Integrate(double x0, double dx, const vector<double> &ys);
void fit_gumbel(double x0, double dx, const vector<double> &ys,
  double &Mu, double &Beta);

static void MakeBins(const vector<float> &TSs, 
  vector<float> &logTSs, vector<uint> &Bins,
  double &x0, double &dx)
	{
	const uint N = SIZE(TSs);

// Sort to facilitate deleting outliers
	vector<float> SortedTSs = TSs;
	sort(SortedTSs.begin(), SortedTSs.end());
	for (uint j = NOUTLIERS; j + NOUTLIERS < N; ++j)
		{
		float TS = SortedTSs[j];
		if (TS > 0)
			{
			float logTs = -logf(TS);
			logTSs.push_back(logTs);
			}
		}
	Binner<float> B(logTSs, NBINS, 0);
	Bins = B.GetBins();

	x0 = B.GetBinMid(0);
	double Mid0 = B.GetBinMid(0);
	double Mid1 = B.GetBinMid(1);
	dx = Mid1 - Mid0;
	}

static void Makexys(const vector<uint> &Bins, double x0, double dx,
  vector<double> &ys)
	{
	ys.clear();
	asserta(SIZE(Bins) == NBINS);
	float Sumy = 0;

	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		{
		uint32_t n = Bins[Bin];
	// trim outliers near x=0, should not be curve-fitted
		if (Bin < 3)
			n = 0;
		float unnormalized_y = float(n);
		ys.push_back(unnormalized_y);
		Sumy += unnormalized_y;
		}

// Normalize so that integral over PDF is 1
	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		ys[Bin] /= (Sumy*dx);

	double S = Integrate(x0, dx, ys);
	asserta(S >= 0.99 && S <= 1.01);
	}

void cmd_calibrate5()
	{
	Die("Not supported");
#if 0
	optset_sensitive = true;
	opt_sensitive = true;
	asserta(!optset_fast);

	if (optset_output)
		Die("Use -calib_output");

	string &QCalFN = g_Arg1;
	if (g_Arg1 == ".")
#ifdef _MSC_VER
		QCalFN = "c:/src/reseek_scop40/reseek_db/scop40_family.cal";
#else
		QCalFN = "/c/src/reseek_scop40/reseek_db/scop40_family.cal";
#endif
	else
		QCalFN = g_Arg1;

	const string &DBFN = opt_db;
	DBSearcher DBS;
	DBS.LoadDB(DBFN);
	const uint NQ = DBS.GetDBChainCount();
	ProgressLog("NQ=%u\n", NQ);
	set<string> Labels;
	for (uint i = 0; i < NQ; ++i)
		{
		const PDBChain &Chain = *DBS.m_DBChains[i];
		const string &Label = Chain.m_Label;
		if (Labels.find(Label) != Labels.end())
			Die("Dupe label >%s", Label.c_str());
		Labels.insert(Label);
		}

	DSSParams Params;
	Params.SetFromCmdLine(DBS.GetDBSize());

	DBS.m_CollectTestStats = true;
	DBS.Setup();
	DBS.RunSelf();
	DSSAligner::Stats();

	FILE *fOut = CreateStdioFile(opt_calib_output);
	const uint DBChainCount = DBS.GetDBChainCount();
	for (uint Idx = 0; Idx < DBChainCount; ++Idx)
		{
		ProgressStep(Idx, DBChainCount, "Calibrating");

		asserta(Idx < SIZE(DBS.m_TestStatsVec));
		const vector<float> &TSs = DBS.m_TestStatsVec[Idx];
		vector<float> logTSs;
		vector<uint> Bins;
		double x0, dx;
		MakeBins(TSs, logTSs, Bins, x0, dx);
		vector<double> ys;
		Makexys(Bins, x0, dx, ys);
		double Mu, Beta;
		fit_gumbel(x0, dx, ys, Mu, Beta);
		if (fOut != 0)
			{
			const PDBChain &Chain = *DBS.m_DBChains[Idx];
			const char *Label = Chain.m_Label.c_str();
			fprintf(fOut, "%s\t%.3g\t%.3g\n", Label, Mu, Beta);
			}
		}
	CloseStdioFile(fOut);
#endif
	}
