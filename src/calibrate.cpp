#include "myutils.h"
#include "calibratesearcher.h"
#include "mx.h"
#include "pdbchain.h"
#include "xdpmem.h"
#include "cigar.h"
#include "binner.h"
#include "fast_erf.h"
#include <thread>
#include <cmath>
#include "timing.h"

double normal(double mu, double sigma, double x)
	{
	static double PI = 3.14159265358979323846;
	double two_sigma2 = 2*sigma*sigma;
	double bottom = sqrt(2*PI*two_sigma2);
	double x_minus_mu2 = (x - mu)*(x - mu);
	double y = (1/bottom)*exp(-x_minus_mu2/two_sigma2);
	return y;
	}

// https://en.wikipedia.org/wiki/Q-function
// Q(x) is the probability that a standard normal variable
//   takes a value larger than x.
double Q_func_standard(double x)
	{
	static const double SQRT2 = sqrt(2);
	//double Q = 0.5 - 0.5*erf(x/SQRT2);
	//double Q = 0.5 - 0.5*SUN_erf(x/SQRT2);
	double Q = 0.5 - 0.5*fast_erf(x/SQRT2);
	return Q;
	}

double Q_func(double x, double mu, double sigma)
	{
	double x_standard = (mu - x)/sigma;
	return Q_func_standard(x_standard);
	}

float GetBinSize(float MinValue, float MaxValue, uint BinCount)
	{
	asserta(BinCount > 1);
	float Range = MaxValue - MinValue;
	asserta(Range > 0);
	float Size = Range/(BinCount - 1);
	return Size;
	}

uint ValueToBin(float MinValue, float MaxValue,
 uint BinCount, float Value)
	{
	if (Value > MaxValue)
		Value = MaxValue;
	if (Value < MinValue)
		Value = MinValue;
	float Range = MaxValue - MinValue;
	asserta(Range > 0);
	float r = (MaxValue - Value)/Range;
	asserta(r >= 0 && r <= 1);
	uint Bin = uint(r*(BinCount-1));
	asserta(Bin < BinCount);
	return Bin;
	}

void CalibrateSearcher::OnSetup()
	{
	m_TestStatsVec.clear();
	m_TestStatsVec.resize(m_ChainCount);
	}

void CalibrateSearcher::OnAln(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA)
	{
	asserta(ChainIndex1 < SIZE(m_TestStatsVec));
	vector<float> &v = m_TestStatsVec[ChainIndex1];
	v.push_back(DA.m_TestStatisticAB);
	}

void CalibrateSearcher::SetAllAccum()
	{
	m_AllAccum.clear();
	m_AllAccum.resize(NBINS, UINT_MAX);
	asserta(SIZE(m_AllBins) == NBINS);
	uint32_t Sum = 0;
	for (uint Bin = 0; Bin < NBINS; ++Bin)
		{
		Sum += m_AllBins[Bin];
		m_AllAccum[Bin] = Sum;
		}
	}

void CalibrateSearcher::SetAllBins()
	{
	asserta(m_ptrAllBinner == 0);
	const vector<vector<float> > &TSV = m_TestStatsVec;
	vector<float> TSs;
	float MaxTS = 0;
	const uint N = SIZE(TSV);
	for (uint i = 0; i < N; ++i)
		{
		vector<float> TSVi = TSV[i];
		sort(TSVi.begin(), TSVi.end());
		const uint n = SIZE(TSVi);
		for (uint j = NOUTLIERS; j < n; ++j)
			{
			float TS = TSVi[j];
			if (TS > 0)
				{
				float logTs = -logf(TS);
				MaxTS = max(logTs, MaxTS);
				TSs.push_back(logTs);
				}
			}
		}
	m_ptrAllBinner = new Binner<float>(TSs, NBINS, 0);
	m_AllBins = m_ptrAllBinner->GetBins();
	asserta(SIZE(m_AllBins) == NBINS);
	}

void CalibrateSearcher::FitNormal()
	{
	asserta(SIZE(m_AllBins) == NBINS);
	float SumValue = 0;
	uint32_t Sumn = 0;
	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		{
		uint32_t n = m_ptrAllBinner->GetCount(Bin);
		float Mid = m_ptrAllBinner->GetBinMid(Bin);
		SumValue += n*Mid;
		Sumn += n;
		}
	m_NormalMeanMinusLogTS = SumValue/Sumn;

	float Sum2 = 0;
	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		{
		uint32_t n = m_ptrAllBinner->GetCount(Bin);
		float Mid = m_ptrAllBinner->GetBinMid(Bin);
		float Diff = Mid - m_NormalMeanMinusLogTS;
		Sum2 += n*Diff*Diff;
		}
	m_NormalSigmaMinusLogTS = sqrtf(Sum2/Sumn);
	}

void CalibrateSearcher::Calibrate(uint ChainIndex, float &Mean, float &Sigma)
	{
	Mean = FLT_MAX;
	Sigma = FLT_MAX;
	asserta(ChainIndex < SIZE(m_TestStatsVec));
	 vector<float> TSVi = m_TestStatsVec[ChainIndex];
	float MaxTS = 0;
	sort(TSVi.begin(), TSVi.end());
	const uint n = SIZE(TSVi);
	vector<float> TSs;
	for (uint j = NOUTLIERS; j < n; ++j)
		{
		float TS = TSVi[j];
		if (TS > 0)
			{
			float logTs = -logf(TS);
			MaxTS = max(logTs, MaxTS);
			TSs.push_back(logTs);
			}
		}
	Binner<float> B(TSs, NBINS, 0);
	const vector<uint32_t> &Bins = B.GetBins();

	asserta(SIZE(Bins) == NBINS);
	float SumValue = 0;
	uint32_t Sumn = 0;
	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		{
		uint32_t n = B.GetCount(Bin);
		float Mid = B.GetBinMid(Bin);
		SumValue += n*Mid;
		Sumn += n;
		}
	Mean = SumValue/Sumn;

	float Sum2 = 0;
	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		{
		uint32_t n = B.GetCount(Bin);
		float Mid = B.GetBinMid(Bin);
		float Diff = Mid - Mean;
		Sum2 += n*Diff*Diff;
		}
	Sigma = sqrtf(Sum2/Sumn);
	}

void CalibrateSearcher::ScanAll()
	{
	m_MaxTS = 0;

// All TestStats after discarding top NOUTLIERS scores for each chain.
	m_AllTSs.clear();

	const vector<vector<float> > &TSV = m_TestStatsVec;
	const uint N = SIZE(TSV);
	for (uint i = 0; i < N; ++i)
		{
		vector<float> TSVi = TSV[i];
		sort(TSVi.begin(), TSVi.end());
		const uint n = SIZE(TSVi);
		for (uint j = NOUTLIERS; j < n; ++j)
			{
			float TS = TSVi[j];
			if (TS > 0)
				{
				m_MaxTS = max(TS, m_MaxTS);
				m_AllTSs.push_back(TS);
				}
			}
		}
	}

void CalibrateSearcher::WriteBins(FILE *f) const
	{
	if (f == 0)
		return;
	uint32_t BinCount = m_ptrAllBinner->GetBinCount();
	fprintf(f, "TS\tMid\tN\tAN\tP\tFit\n");
	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		{
		uint32_t n = m_ptrAllBinner->GetCount(Bin);
		asserta(m_AllBins[Bin] == n);
		uint32_t an = m_AllAccum[Bin];
		float Mid = m_ptrAllBinner->GetBinMid(Bin);
		double TS = exp(-Mid);
		double Fit = normal(m_NormalMeanMinusLogTS, m_NormalSigmaMinusLogTS, Mid);
		double P = Q_func(Mid, m_NormalMeanMinusLogTS, m_NormalSigmaMinusLogTS);
		fprintf(f, "%.3g\t%.3f\t%u\t%u\t%.3g\t%.3g\n", TS, Mid, n, an, P, Fit);
		}
	}

void cmd_calibrate()
	{
	const string &QCalFN = g_Arg1;
	const string &DBFN = opt_db;
	CalibrateSearcher DBS;
	DSSParams Params;
	Params.SetFromCmdLine();
	DBS.ReadChains(QCalFN, DBFN);
	Params.m_DBSize = (float) DBS.GetDBSize();
	if (optset_dbsize)
		Params.m_DBSize = (float) opt_dbsize;
	FILE *fOut = CreateStdioFile(opt_output);

	DBS.Setup(Params);
	DBS.Run();
	DSSAligner::Stats();
	DBS.ScanAll();
	DBS.SetAllBins();
	DBS.SetAllAccum();
	DBS.FitNormal();
	DBS.WriteBins(fOut);
	CloseStdioFile(fOut);

	Log("Mean %.3g, stddev %.3g\n", DBS.m_NormalMeanMinusLogTS, DBS.m_NormalSigmaMinusLogTS);
	}
