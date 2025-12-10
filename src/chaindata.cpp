#include "myutils.h"
#include "chaindata.h"

const float *ChainData::GetDistMxPtr()
	{
	if (m_DistMxPtr.m_Filled)
		return m_DistMxPtr.m_Data;

	float *DistMxPtr = m_DistMxPtr.Alloc(m_L*m_L);
	for (uint i = 0; i < m_L; ++i)
		{
		DistMxPtr[i*m_L + i] = 0;
		for (uint j = i + 1; j < m_L; ++j)
			{
			float d = m_Chain->GetDist(i, j);
			DistMxPtr[i*m_L + j] = d;
			DistMxPtr[j*m_L + i] = d;
			}
		}
	m_DistMxPtr.Filled();
	return DistMxPtr;
	}
