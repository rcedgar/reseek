#pragma once

#include "dbsearcher.h"

class DistMxSearcher : public DBSearcher
	{
public:
	float m_MaxTS = 0;
	FILE *m_fDistMx = 0;
	uint m_ChainCount = 0;

public:
	virtual void OnSetup();
	virtual void OnAln(DSSAligner &DA, bool Up);
	};
