#include "myutils.h"
#include "calibratesearcher3.h"
#include "binner.h"

void CalibrateSearcher3::OnSetup()
	{
	m_TestStatsVec.clear();
	m_TestStatsVec.resize(m_ChainCount);
	}

void CalibrateSearcher3::OnAln(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA)
	{
	asserta(ChainIndex1 < SIZE(m_TestStatsVec));
	vector<float> &v = m_TestStatsVec[ChainIndex1];
	v.push_back(DA.m_TestStatisticA);
	}

void CalibrateSearcher3::OnAlnBA(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA)
	{
	asserta(ChainIndex2 < SIZE(m_TestStatsVec));
	vector<float> &v = m_TestStatsVec[ChainIndex2];
	v.push_back(DA.m_TestStatisticB);
	}
