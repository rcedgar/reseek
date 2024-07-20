#include "myutils.h"
#include "calibratesearcher.h"
#include "chainreader2.h"
#include "binner.h"
#include "timing.h"
#include "sort.h"
#include <set>

/***
* -calibrate3
* Calculate per-chain count bins for TS values.
***/

void CalibrateSearcher::WriteTopFPsBottomTPs(FILE *f) const
	{
	if (f == 0)
		return;

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
		if (NTP > 0)
			QuickSortInPlaceDesc(TPs.data(), SIZE(TPs));
		if (NFP > 0)
			QuickSortInPlaceDesc(FPs.data(), SIZE(FPs));
		fprintf(f, "%s\ttps=", Dom.c_str());
		for (uint i = 0; i < min(NFP, 10u); ++i)
			{
			if (i > 0)
				fprintf(f, ",");
			fprintf(f, "%.3g", TPs[i]);
			}
		fprintf(f, "\tfps=");
		for (uint i = 0; i < min(NFP, 10u); ++i)
			{
			if (i > 0)
				fprintf(f, ",");
			fprintf(f, "%.3g", FPs[i]);
			}
		fprintf(f, "\n");
		}
	}

 //32 bins, TS range 0.0 to 0.3, written to -calib_output calibrate3.tsv
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
