#include "myutils.h"
#include "calibratesearcher3.h"
#include "binner.h"
#include "timing.h"
#include <set>

/***
* Calibrate E-value by measuring score 
* distribution (unsupervised).
***/
void cmd_calibrate3()
	{
	asserta(optset_output);
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
	CalibrateSearcher3 DBS;
	DBS.ReadChains(QCalFN, DBFN);
	const uint NQ = DBS.GetQueryCount();
	ProgressLog("NQ=%u\n", NQ);
	set<string> Labels;
	for (uint i = 0; i < NQ; ++i)
		{
		const PDBChain &Chain = *DBS.m_Chains[i];
		const string &Label = Chain.m_Label;
		if (Labels.find(Label) != Labels.end())
			Die("Dupe label >%s", Label.c_str());
		Labels.insert(Label);
		}

	DSSParams Params;
	Params.SetFromCmdLine(DBS.GetDBSize());

	const uint SAMPLE = 32;
	const uint NBINS = 32;
	DBS.Setup(Params);
	DBS.Run();
	DSSAligner::Stats();

	FILE *fOut = CreateStdioFile(opt_output);
	vector<float> Cutoffs;
	float SumCutoffs = 0;
	for (uint i = 0; i < NQ; ++i)
		{
		ProgressStep(i, NQ, "Cutoffs\n");
		const PDBChain &Chain = *DBS.m_Chains[i];
		const vector<float> &TSs_i = DBS.m_TestStatsVec[i];
		Binner<float> B(TSs_i, NBINS, -0.1f, 0.3f);
		vector<uint> Bins;
		B.GetAccumBinsReverse(Bins);
		float Cutoff = B.GetCutoff_Fract(0.5);
		Cutoffs.push_back(Cutoff);
		SumCutoffs += Cutoff;
		const char *Label = Chain.m_Label.c_str();
		fprintf(fOut, "%s\t%.3g\n", Label, Cutoff);
		}
	CloseStdioFile(fOut);
	fOut = 0;
	float MeanCutoff = SumCutoffs/NQ;
	ProgressLog("MeanCutoff %.3g\n", MeanCutoff);

//////////////////////////////////////////////////////
//	Write small sample of per-query distributions to 
//   log file for plotting
//////////////////////////////////////////////////////

// Mid-point TS values for bins
	vector<float> Mids;
	const vector<float> &TSs_0 = DBS.m_TestStatsVec[0];
	Binner<float> B(TSs_0, NBINS, -0.1f, 0.3f);
	vector<uint> Bins;
	B.GetAccumBinsReverse(Bins);
	for (uint Bin = 0; Bin < NBINS; ++Bin)
		Mids.push_back(B.GetBinMid(Bin));
 
	vector<vector<uint> > AccumBinsVec(SAMPLE);
	for (uint i = 0; i < SAMPLE; ++i)
		{
		const vector<float> &TSs_i = DBS.m_TestStatsVec[i];
		Binner<float> B(TSs_i, NBINS, -0.1f, 0.3f);
		B.GetAccumBinsReverse(AccumBinsVec[i]);
		}

	for (uint Bin = 0; Bin < NBINS; ++Bin)
		{
		fprintf(g_fLog, "%.4g", Mids[Bin]);
		for (uint i = 0; i < SAMPLE; ++i)
			{
			uint n = AccumBinsVec[i][Bin];
			if (n == 0)
				fprintf(g_fLog, "\t");
			else
				fprintf(g_fLog, "\t%u", n);
			}
		fprintf(g_fLog, "\n");
		}

	fprintf(g_fLog, "Cutoffs");
	for (uint i = 0; i < SAMPLE; ++i)
		fprintf(g_fLog, "\t%.3g", Cutoffs[i]);
	fprintf(g_fLog, "\n");
	fprintf(g_fLog, "MeanC\t%.3g\n", MeanCutoff);
//////////////////////////////////////////////////////
//  End sample
//////////////////////////////////////////////////////
	}
