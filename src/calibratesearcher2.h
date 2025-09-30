#if 0 // @@DELETE
#pragma once

#include "scop40bench.h"
#include "binner.h"

/***
* Calibrate distribution of FP rates as function of test statistic
* on SCOP40 or similar.
***/
class CalibrateSearcher2 : public SCOP40Bench
	{
public:
	Binner<float> *m_ptrAllBinner = 0;  // binned by -log(TS)
	float m_MaxTS = 0;
	vector<float> m_AllTSs; // outliers discarded
	vector<uint32_t> m_AllBins;
	vector<uint32_t> m_AllAccum;
	};
#endif