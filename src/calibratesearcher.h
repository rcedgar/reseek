#pragma once

#include "scop40bench.h"
#include "binner.h"

/***
* Calibrate distribution of FP errors on one_per_sf or one_per_fold
* Gumbel fits well for the bulk of the distribution, but the fit is
* poor in both tails which causes E-value estimates to diverge in
* the most relevant ranges for practice.
***/
class CalibrateSearcher : public SCOP40Bench
	{
public:
	Binner<float> *m_ptrAllBinner = 0;  // binned by -log(TS)
	float m_MaxTS = 0;
	vector<float> m_AllTSs; // outliers discarded
	vector<uint32_t> m_AllBins;
	vector<uint32_t> m_AllAccum;
	vector<vector<float> > m_TestStatsVec;

	double m_NormalMean = 0;
	double m_NormalSigma = 0;

	double m_GumbelMu = 0;
	double m_GumbelBeta = 0;

// For fitting to distribution
//   x values are -log(TS)
//   y values are bin frequencies (sum(y) = 1)
	double m_x0 = 0;
	double m_dx = 0;
	vector<double> m_ys;

public:
	virtual void OnSetup();
	virtual void OnAln(DSSAligner &DA, bool Up);

public:
// All-domain calibration
	void ScanAll();
	void SetAllBins();
	void Setxys();
	void FitNormal();
	void FitGumbel();
	void SetAllAccum();
	void WriteBins(FILE *f) const;

// Per-domain calibration
	void WriteSlopeCalibOutput(FILE *f,
	  uint BinCount, float TSlo, float TShi) const;
	void WriteTopFPsBottomTPs(FILE *f, uint n) const;
	};

static const uint NOUTLIERS = 3;
static const uint NBINS = 101;
