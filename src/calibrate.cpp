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

double gumbel(double mu, double beta, double x);
double gumbel_cdf(double mu, double beta, double x);

double Integrate(double x0, double dx,
  const vector<double> &ys)
	{
	const uint N = SIZE(ys);
	double Sum = 0;
	double x = x0;
	for (uint i = 0; i < N; ++i)
		{
		double y = ys[i];
		Sum += y*dx;
		x += dx;
		}
	return Sum;
	}

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

void CalibrateSearcher::Setxys()
	{
	m_ys.clear();
	asserta(SIZE(m_AllBins) == NBINS);
	float Sumy = 0;

	m_x0 = m_ptrAllBinner->GetBinMid(0);
	double Mid0 = m_ptrAllBinner->GetBinMid(0);
	double Mid1 = m_ptrAllBinner->GetBinMid(1);

	m_dx = Mid1 - Mid0;
	double x = m_x0;

	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		{
		uint32_t n = m_ptrAllBinner->GetCount(Bin);
	// trim outliers near x=0, should not be curve-fitted
		if (Bin < 10)
			n = 0;
		float unnormalized_y = float(n);
		m_ys.push_back(unnormalized_y);
		Sumy += unnormalized_y;
		x += m_dx;
		}

// Normalize so that integral over PDF is 1
	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		m_ys[Bin] /= (Sumy*m_dx);

	double S = Integrate(m_x0, m_dx, m_ys);
	asserta(S >= 0.99 && S <= 1.01);
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

void CalibrateSearcher::FitGumbel()
	{
	void fit_gumbel(double x0, double dx,
	  const vector<double> &ys,
	  double &Mu, double &Beta);

	fit_gumbel(m_x0, m_dx, m_ys, m_GumbelMu, m_GumbelBeta);

	Log("Gumbel: Mu %.3g, Beta %.3g\n",
	  m_GumbelMu, m_GumbelBeta);
	}

void CalibrateSearcher::FitNormal()
	{
	asserta(SIZE(m_ys) == NBINS);

	double Sumxy = 0;
	double Sumy = 0;
	double x = m_x0;
	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		{
		double y = m_ys[Bin];
		Sumy += y;
		Sumxy += x*y;
		x += m_dx;
		}
	m_NormalMean = Sumxy/Sumy;

	double Sum2 = 0;
	x = m_x0;
	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		{
		double y = m_ys[Bin];
		double Diff = x - m_NormalMean;
		Sum2 += y*Diff*Diff;
		x += m_dx;
		}
	m_NormalSigma = sqrt(Sum2);

	Log("Normal: Mean %.3g, stddev %.3g\n",
	  m_NormalMean, m_NormalSigma);
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

	fprintf(f, "Bin");
	fprintf(f, "\tTS");
	fprintf(f, "\tMid");
	fprintf(f, "\tx");
	fprintf(f, "\tn");
	fprintf(f, "\tan");
	fprintf(f, "\ty");
	fprintf(f, "\ty_fit");
	fprintf(f, "\tx0=%.3g\tdx=%.3g\n", m_x0, m_dx);

	double x = m_x0;
	for (uint32_t Bin = 0; Bin < NBINS; ++Bin)
		{
		uint32_t n = m_ptrAllBinner->GetCount(Bin);
		asserta(m_AllBins[Bin] == n);
		uint32_t an = m_AllAccum[Bin];
		double y = m_ys[Bin];
		float Mid = m_ptrAllBinner->GetBinMid(Bin);
		double TS = exp(-x);
		double Fit = gumbel(m_GumbelMu, m_GumbelBeta, x);
		double P = gumbel_cdf(m_GumbelMu, m_GumbelBeta, x);

		fprintf(f, "%u", Bin);
		fprintf(f, "\t%.3g", TS);
		fprintf(f, "\t%.3g", Mid);
		fprintf(f, "\t%.3g", x);
		fprintf(f, "\t%u", n);
		fprintf(f, "\t%u", an);
		fprintf(f, "\t%.3g", y);
		fprintf(f, "\t%.3g", Fit);
		fprintf(f, "\t%.3g", P);
		fprintf(f, "\n");

		x += m_dx;
		}
	}

static double GetE1(uint DBSize, double TS)
	{
	const float Slope = -6.6f; // -7.3f;
	const float Intercept = 6.1f;
	double logNF = Slope*TS+ Intercept;
	double NF = pow(10, logNF);
	double Evalue = NF*DBSize/1e8f;
	return Evalue;
	}

// Gumbel: Mu 2.5, Beta 0.613
static double GetE2(uint DBSize, double TS)
	{
	double x = -log(TS);
	double P = gumbel_cdf(2.5, 0.613, x);
	double E = P*DBSize;
	return E;
	}

void cmd_evalue_table()
	{
	double TS = 1.0;
	uint DBSize = 10000;
	for (;;)
		{
		if (TS < 0.0001)
			break;
		TS *= 0.9;
		double E1 = GetE1(DBSize, TS);
		double E2 = GetE2(DBSize, TS);
		Log("%.3g\t%.3g\t%.3g\n", TS, E1, E2);
		}
	}

void cmd_calibrate()
	{
	const string &QCalFN = g_Arg1;
	const string &DBFN = opt_db;
	CalibrateSearcher DBS;
	DBS.ReadChains(QCalFN, "");

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
