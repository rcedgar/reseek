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

void DSSAligner::AlignBags(const ChainBag &BagA,
						   const ChainBag &BagB)
	{
	ClearAlign();
	if (DoMKF_Bags(BagA, BagB))
		{
		m_MKF.SetBagQ(BagA);
		m_MKF.AlignBag(BagB);
		PostAlignMKF();
		return;
		}

	float Omega = m_Params->m_Omega;
	if (Omega > 0)
		{
		float MuScore = AlignMuParaBags(BagA, BagB);
		if (MuScore < Omega)
			return;
		}
	SetSMx_NoRev(*m_Params, *BagA.m_ptrProfile, *BagB.m_ptrProfile);
	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	uint Leni, Lenj;
	m_AlnFwdScore = SWFast(m_Mem, GetSMxData(), LA, LB,
	  m_Params->m_GapOpen, m_Params->m_GapExt,
	  m_LoA, m_LoB, Leni, Lenj, m_Path);

	CalcEvalue();
	}
