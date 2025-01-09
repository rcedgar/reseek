#include "myutils.h"
#include "chainbag.h"
#include "dssaligner.h"

bool DSSAligner::DoMKF_Bags(const ChainBag &BagA,
							const ChainBag &BagB) const
	{
	if (BagA.m_ptrMuLetters == 0 || BagB.m_ptrMuLetters == 0)
		return false;
	if (BagA.m_ptrMuLetters == 0 || BagB.m_ptrMuLetters == 0)
		return false;

	uint LA = BagA.m_ptrChain->GetSeqLength();
	uint LB = BagB.m_ptrChain->GetSeqLength();
	if (LA >= m_Params->m_MKFL)
		return true;
	if (LB >= m_Params->m_MKFL)
		return true;
	return false;
	}

void DSSAligner::AlignBags(const ChainBag &BagQ,
						   const ChainBag &BagT)
	{
	ClearAlign();
	if (DoMKF_Bags(BagQ, BagT))
		{
		m_MKF.SetBagQ(BagQ);
		m_MKF.AlignBag(BagT);
		PostAlignMKF();
		return;
		}
	}
