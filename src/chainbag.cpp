#include "myutils.h"
#include "chainbag.h"
#include "dssaligner.h"
#include "parasail.h"

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

	m_ChainA = BagA.m_ptrChain;
	m_ChainB = BagB.m_ptrChain;

	m_ProfileA = BagA.m_ptrProfile;
	m_ProfileB = BagB.m_ptrProfile;

	m_SelfRevScoreA = BagA.m_SelfRevScore;
	m_SelfRevScoreB = BagB.m_SelfRevScore;

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
	const uint LA = BagA.m_ptrChain->GetSeqLength();
	const uint LB = BagB.m_ptrChain->GetSeqLength();

	uint Leni, Lenj;
	m_AlnFwdScore = SWFast(m_Mem, GetSMxData(), LA, LB,
	  m_Params->m_GapOpen, m_Params->m_GapExt,
	  m_LoA, m_LoB, Leni, Lenj, m_Path);

	CalcEvalue();
	}

void ChainBag::Validate(const char *Msg) const
	{
	if (m_ptrChain == 0)
		{
		asserta(m_ptrProfile == 0);
		asserta(m_ptrMuLetters == 0);
		asserta(m_ptrMuKmers == 0);
		asserta(m_ptrProfPara == 0);
		asserta(m_ptrProfParaRev == 0);
		asserta(m_ptrKmerHashTableQ == 0);
		return;
		}
	const uint L = m_ptrChain->GetSeqLength();
	asserta(SIZE(*m_ptrMuLetters) == L);
	asserta(SIZE(*m_ptrMuKmers) + 2 == L);
	const parasail_profile_t *Prof = (const parasail_profile_t *) m_ptrProfPara;
	const parasail_profile_t *ProfRev = (const parasail_profile_t *) m_ptrProfParaRev;
	if (Prof != 0)
		asserta(Prof->s1Len == L);
	if (ProfRev != 0)
		asserta(ProfRev->s1Len == L);
	}
