#pragma once

#include "dbsearcher.h"
#include "binner.h"

class CalibrateSearcher3 : public DBSearcher
	{
public:
	vector<vector<float> > m_TestStatsVec;

public:
	virtual void OnSetup();
	virtual void OnAln(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA);
	virtual void OnAlnBA(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA);
	};

