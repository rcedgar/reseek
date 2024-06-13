#pragma once

#include "scop40bench.h"
#include "binner.h"

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
	virtual void OnAln(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA);
	virtual void OnAlnBA(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA);

public:
	void ScanAll();
	void SetAllBins();
	void Setxys();
	void FitNormal();
	void FitGumbel();
	void SetAllAccum();
	void WriteBins(FILE *f) const;
	};

static const uint NOUTLIERS = 3;
static const uint NBINS = 101;
