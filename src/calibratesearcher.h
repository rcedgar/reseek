#pragma once

#include "dbsearcher.h"

class CalibrateSearcher : public DBSearcher
	{
public:
	vector<vector<float> > m_TestStatsVec;

public:
	virtual void OnSetup();
	virtual void OnAln(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA);
	};
