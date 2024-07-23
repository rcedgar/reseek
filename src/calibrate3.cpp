#include "myutils.h"
#include "calibratesearcher.h"
#include "chainreader2.h"
#include "binner.h"
#include "timing.h"
#include "sort.h"
#include <set>

/***
* -calibrate3
* Calculate per-chain TP and FP count bins for TS values.
***/
void CalibrateSearcher::WriteTP_FP_TSBins(FILE *f, uint BinCount,
  float TSlo, float TShi) const
	{
	if (f == 0)
		return;

// Mid-point TS values for bins
	vector<float> Mids;
	const vector<float> &TSs_0 = m_TestStatsVec[0];
	Binner<float> B0(TSs_0, BinCount, TSlo, TShi);
	const vector<uint> &Bins0 = B0.GetBins();
	fprintf(f, "TS");
	for (uint Bin = 0; Bin < BinCount; ++Bin)
		fprintf(f, "\tTP%.3g", B0.GetBinMid(Bin));
	for (uint Bin = 0; Bin < BinCount; ++Bin)
		fprintf(f, "\tFP%.3g", B0.GetBinMid(Bin));
	fprintf(f, "\n");

	const uint DomCount = GetDomCount();
	vector<vector<float> > DomIdxToTPs(DomCount);
	vector<vector<float> > DomIdxToFPs(DomCount);
	const uint HitCount = GetHitCount();
	for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
		{
		uint DomIdx1 = m_DomIdx1s[HitIdx];
		uint DomIdx2 = m_DomIdx2s[HitIdx];
		if (DomIdx1 == DomIdx2)
			continue;
		float TS = m_TSs[HitIdx];
		bool T = IsT(DomIdx1, DomIdx2);
		if (T)
			DomIdxToTPs[DomIdx1].push_back(TS);
		else
			DomIdxToFPs[DomIdx1].push_back(TS);
		}

	for (uint DomIdx = 0; DomIdx < DomCount; ++DomIdx)
		{
		const string &Dom = m_Doms[DomIdx];
		vector<float> &TPs = DomIdxToTPs[DomIdx];
		vector<float> &FPs = DomIdxToFPs[DomIdx];
		const uint NTP = SIZE(TPs);
		const uint NFP = SIZE(FPs);
		Binner<float> TPB(TPs, BinCount, TSlo, TShi);
		Binner<float> FPB(FPs, BinCount, TSlo, TShi);
		const vector<uint> &TPBins = TPB.GetBins();
		const vector<uint> &FPBins = FPB.GetBins();
		fprintf(f, "%s", Dom.c_str());
		for (uint Bin = 0; Bin < BinCount; ++Bin)
			{
			uint n = TPBins[Bin];
			fprintf(f, "\t%u", n);
			}
		for (uint Bin = 0; Bin < BinCount; ++Bin)
			{
			uint n = FPBins[Bin];
			fprintf(f, "\t%u", n);
			}
		fprintf(f, "\n");
		}
	}

void CalibrateSearcher::WriteTSBins(FILE *f, uint BinCount,
  float TSlo, float TShi) const
	{
	if (f == 0)
		return;

// Mid-point TS values for bins
	vector<float> Mids;
	const vector<float> &TSs_0 = m_TestStatsVec[0];
	Binner<float> B(TSs_0, BinCount, TSlo, TShi);
	const vector<uint> &Bins = B.GetBins();
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

	const uint BIN_COUNT = 20;
	const float MIN_TS = 0.0f;
	const float MAX_TS = 0.8f;

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
	DBS.RunSelf();

	DSSAligner::Stats();

	if (optset_calib_output)
		{
		FILE *fOut = CreateStdioFile(opt_calib_output);
		DBS.WriteTSBins(fOut, BIN_COUNT, MIN_TS, MAX_TS);
		CloseStdioFile(fOut);
		}

	if (optset_calib_output2)
		{
		FILE *fOut = CreateStdioFile(opt_calib_output2);
		DBS.m_Level = string("sf");
		DBS.WriteTP_FP_TSBins(fOut, BIN_COUNT, MIN_TS, MAX_TS);
		CloseStdioFile(fOut);
		}
	}
