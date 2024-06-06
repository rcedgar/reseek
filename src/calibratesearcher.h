#pragma once

#include "dbsearcher.h"
#include "binner.h"

class CalibrateSearcher : public DBSearcher
	{
public:
	Binner<float> *m_ptrAllBinner = 0;
	float m_MaxTS = 0;
	vector<float> m_AllTSs; // outliers discarded
	vector<uint32_t> m_AllBins;
	vector<vector<float> > m_TestStatsVec;
	float m_AllMean = 0;
	float m_AllSigma = 0;

public:
	virtual void OnSetup();
	virtual void OnAln(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA);

public:
	void ScanAll();
	void SetAllBins();
	void SetAllMeanStdDev();
	void LogAllBins() const;
	void Calibrate(uint ChainIndex, float &Mean, float &Sigma);
	};

static const uint NOUTLIERS = 3;
static const uint NBINS = 101;
