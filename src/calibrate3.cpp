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
void CalibrateSearcher::WriteTP_FP_TSBins(const string &FN,
  uint BinCount, float TSlo, float TShi) const
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);

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
	CloseStdioFile(f);
	}

void CalibrateSearcher::GetTSBins(uint BinCount, float TSlo, float TShi,
  vector<float> &Mids, vector<string> &Labels, vector<vector<uint> > &BinsVec) const
	{
// Mid-point TS values for bins
	Mids.clear();
	const vector<float> &TSs_0 = m_TestStatsVec[0];
	Binner<float> B(TSs_0, BinCount, TSlo, TShi);
	for (uint Bin = 0; Bin < BinCount; ++Bin)
		Mids.push_back(B.GetBinMid(Bin));

// Vectors indexed by DB chain idx
	Labels.clear();
	BinsVec.clear();
	const uint NR = GetDBChainCount();
	BinsVec.resize(NR);
	for (uint i = 0; i < NR; ++i)
		{
		const string &Label = m_DBChains[i]->m_Label;
		Labels.push_back(Label);

		const vector<float> &TSs_i = m_TestStatsVec[i];
		Binner<float> B(TSs_i, BinCount, TSlo, TShi);
		BinsVec[i] = B.GetBins();
		}
	}

void CalibrateSearcher::WriteTSVec(const string &FN) const
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);
	const uint NR = GetDBChainCount();
	for (uint i = 0; i < NR; ++i)
		{
		const char *Label = m_DBChains[i]->m_Label.c_str();
		const vector<float> &TSs_i = m_TestStatsVec[i];
		const uint N = SIZE(TSs_i);
		fprintf(f, "%s", Label);
		for (uint j = 0; j < N; ++j)
			fprintf(f, "\t%.3g", TSs_i[j]);
		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	}

void CalibrateSearcher::WriteTSBins(const string &FN,
  uint BinCount, float TSlo, float TShi) const
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);

// Mid-point TS values for bins
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
	CloseStdioFile(f);
	}

static uint GetTSBin(float TS)
	{
	TS *= 10;
	if (TS >= 10)
		return 10;
	if (TS <= 0)
		return 0;
	uint Bin = uint(TS + 0.5);
	asserta(Bin <= 10);
	return Bin;
	}

static uint AddNoise1(const uint x, uint Pct)
	{
	uint Counter = 0;
	asserta(Pct >= 0 && Pct <= 100);
	int sign = (randu32()%2 == 0 ? 1 : -1);
	uint Pctd = randu32()%(Pct + 1);
	double d = Pctd/100.0;
	uint noisy_x = uint(x*(1 + sign*d));
	return noisy_x;
	}

static void AddNoise(const vector<uint> &x, 
  uint Pct, vector<uint> &y)
	{
	y.clear();
	const uint N = SIZE(x);
	for (uint i = 0; i < N; ++i)
		y.push_back(AddNoise1(x[i], Pct));
	}

void GetFeatures(const vector<uint> &Bins, vector<float> &Features)
	{
	Features.clear();
	const uint N = SIZE(Bins);
	uint maxn = 0;
	for (uint i = 0; i < N; ++i)
		maxn = max(maxn, Bins[i]);
	float FM = float(maxn) + 1;
	for (uint i = 0; i < N; ++i)
		Features.push_back(Bins[i]/FM);
	}

static void Output(FILE *f,
  const CalibrateSearcher &DBS,
  uint BIN_COUNT, uint NOISE_PCT,
  const vector<vector<uint> > &BinsVec,
  const vector<uint> &XPs, bool T)
	{
	const char Tc = (T ? '1' : '0');
	const uint NXP = SIZE(XPs);
	for (uint i = 0; i < NXP; ++i)
		{
		uint HitIdx = XPs[i];
		uint DomIdx1 = DBS.m_DomIdx1s[HitIdx];
		uint DomIdx2 = DBS.m_DomIdx2s[HitIdx];
		const float TS = DBS.m_TSs[HitIdx];
		int iT = DBS.IsT(DomIdx1, DomIdx2);
		if (T)
			asserta(iT == 1);
		else
			asserta(iT == 0);
		uint ChainIdx2 = DBS.m_HitChainIdxs[HitIdx];
		const string &Label = DBS.m_DBChains[ChainIdx2]->m_Label;

		const vector<uint> &Bins = BinsVec[ChainIdx2];
		asserta(SIZE(Bins) == BIN_COUNT);
		vector<uint> NoisyBins;

		AddNoise(Bins, NOISE_PCT, NoisyBins);
		asserta(SIZE(NoisyBins) == BIN_COUNT);

		vector<float> Features;
		GetFeatures(NoisyBins, Features);
		asserta(SIZE(Features) == BIN_COUNT);

		fprintf(f, "%s", Label.c_str());
		fprintf(f, "\t%c", Tc);
		fprintf(f, "\t%.3g", TS);
		for (uint Bin = 0; Bin < BIN_COUNT; ++Bin)
			fprintf(f, "\t%.3g", Features[Bin]);
		fprintf(f, "\n");
		}
	}

void cmd_calibrate3()
	{
	if (optset_output)
		Die("Use -calib_output");

	const string &DBFN = g_Arg1;

	optset_sensitive = true;
	opt(sensitive) = true;

	optset_minchainlength = true;
	opt(minchainlength) = 1;

	optset_evalue = true;
	opt(evalue) = 10;

	const uint BIN_COUNT = 16;
	const float MIN_TS = 0.0f;
	const float MAX_TS = 1.0f;
	const uint MAX_HITS_PER_TS_BIN_TP = 1000;
	const uint MAX_HITS_PER_TS_BIN_FP = 5000;
	const uint MIN_HITS_PER_TS_BIN = 100;
	const uint NOISE_PCT = 20;

	CalibrateSearcher DBS;
	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast, SCOP40_DBSIZE);
	DBS.m_Params = &Params;

	DBS.LoadDB(DBFN);
	DBS.Setup();
	DBS.RunSelf();
	DSSAligner::Stats();

	DBS.WriteTSBins(opt(calib_output), BIN_COUNT, MIN_TS, MAX_TS);
	DBS.WriteTP_FP_TSBins(opt(calib_output2), BIN_COUNT, MIN_TS, MAX_TS);

	if (!optset_calib_output3)
		return;

	vector<vector<uint> > BinsVec;
	vector<float> Mids;
	vector<string> Labels;
	DBS.GetTSBins(BIN_COUNT, MIN_TS, MAX_TS, Mids, Labels, BinsVec);

	vector<vector<uint> > TSBinToHitIdxs_TP(11);
	vector<vector<uint> > TSBinToHitIdxs_FP(11);
	const uint HitCount = DBS.GetHitCount();
	for (uint HitIdx = 0; HitIdx < HitCount; ++HitIdx)
		{
		uint DomIdx1 = DBS.m_DomIdx1s[HitIdx];
		uint DomIdx2 = DBS.m_DomIdx2s[HitIdx];
		float TS = DBS.m_TSs[HitIdx];
		uint TSBin = GetTSBin(TS);
		if (DBS.IsT(DomIdx1, DomIdx2))
			TSBinToHitIdxs_TP[TSBin].push_back(HitIdx);
		else
			TSBinToHitIdxs_FP[TSBin].push_back(HitIdx);
		}

	for (uint TSBin = 0; TSBin <= 10; ++TSBin)
		{
		vector<uint> &TPs = TSBinToHitIdxs_TP[TSBin];
		vector<uint> &FPs = TSBinToHitIdxs_FP[TSBin];
		const uint NTP = SIZE(TPs);
		const uint NFP = SIZE(FPs);
		const uint N = NTP + NFP;
		asserta(N > 0);
		if (NTP > 0)
			Shuffle(TPs);
		if (NFP > 0)
			Shuffle(FPs);

		double FractionTrue = double(NTP)/double(N);
		uint SampledNTP = uint(FractionTrue*MAX_HITS_PER_TS_BIN_TP);
		uint SampledNFP = uint((1 - FractionTrue)*MAX_HITS_PER_TS_BIN_FP);

		if (SampledNTP < MIN_HITS_PER_TS_BIN)
			SampledNTP = MIN_HITS_PER_TS_BIN;
		if (SampledNFP < MIN_HITS_PER_TS_BIN)
			SampledNFP = MIN_HITS_PER_TS_BIN;

		if (SampledNFP > NFP)
			SampledNFP = NFP;
		if (SampledNTP > NTP)
			SampledNTP = NTP;

		FPs.resize(SampledNFP);
		TPs.resize(SampledNTP);
		}

	for (uint TSBin = 0; TSBin <= 10; ++TSBin)
		{
		const vector<uint> &TPs = TSBinToHitIdxs_TP[TSBin];
		const vector<uint> &FPs = TSBinToHitIdxs_FP[TSBin];
		const uint NTP = SIZE(TPs);
		const uint NFP = SIZE(FPs);
		Log("TSBin %2u  TPs %u FPs %u\n", TSBin, NTP, NFP);
		}

	FILE *f = CreateStdioFile(opt(calib_output3));
	fprintf(f, "Label");
	fprintf(f, "\tT");
	fprintf(f, "\tTS");
	for (uint Bin = 0; Bin < BIN_COUNT; ++Bin)
		fprintf(f, "\tF%u", Bin);
	fprintf(f, "\n");

	for (uint TSBin = 0; TSBin <= 10; ++TSBin)
		{
		const vector<uint> &TPs = TSBinToHitIdxs_TP[TSBin];
		const vector<uint> &FPs = TSBinToHitIdxs_FP[TSBin];
		Output(f, DBS, BIN_COUNT, NOISE_PCT, BinsVec, TPs, true);
		Output(f, DBS, BIN_COUNT, NOISE_PCT, BinsVec, FPs, false);
		}

	CloseStdioFile(f);
	}
