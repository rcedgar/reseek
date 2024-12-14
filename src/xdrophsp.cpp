#include "myutils.h"
#include "dssaligner.h"

float DSSAligner::SubstScore(uint PosA, uint PosB)
	{
	const DSSParams &Params = *m_Params;
	const vector<vector<byte> > &ProfileA = *m_ProfileA;
	const vector<vector<byte> > &ProfileB = *m_ProfileB;
	const uint FeatureCount = Params.GetFeatureCount();
	assert(SIZE(ProfileA) == FeatureCount);
	assert(SIZE(ProfileB) == FeatureCount);
	float Total = 0;
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = m_Params->m_Features[FeatureIdx];
		uint AlphaSize = g_AlphaSizes2[F];
		float **ScoreMx = m_Params->m_ScoreMxs[F];
		const vector<byte> &ProfRowA = ProfileA[FeatureIdx];
		const vector<byte> &ProfRowB = ProfileB[FeatureIdx];
		byte ia = ProfRowA[PosA];
		assert(ia < AlphaSize);
		const float *ScoreMxRow = ScoreMx[ia];

		byte ib = ProfRowB[PosB];
		assert(ib < AlphaSize);
		Total += ScoreMxRow[ib];
		}
	return Total;
	}

float DSSAligner::StaticSubstScore(void *UserData_this, uint PosA, uint PosB)
	{
	DSSAligner *pThis = (DSSAligner *) UserData_this;
	float Score = pThis->SubstScore(PosA, PosB);
	return Score;
	}

float DSSAligner::XDropHSP(uint Loi, uint Loj, uint Len)
	{
	const float Open = m_Params->m_GapOpen;
	const float Ext = m_Params->m_GapExt;
	const float X = float(m_Params->m_MKF_X2);
	//const float X = 8;
	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	uint LoA = Loi + Len/2;
	uint LoB = Loj + Len/2;
	asserta(LoA > 0 && LoB > 0);

	string PathFwd, PathBwd;
	float ScoreFwd = 
		XDropFwd(m_Mem, X, Open, Ext, StaticSubstScore, (void *) this, 
				 LoA, LA, LoB, LB, PathFwd);

	float ScoreBwd = 
		XDropBwd(m_Mem, X, Open, Ext, StaticSubstScore, (void *) this, 
				 LoA-1, LA, LoB-1, LB, PathBwd);
	float TotalScore = ScoreFwd + ScoreBwd;
	if (TotalScore < 10)
		{
		m_XDropPath.clear();
		return 0;
		}

	uint MergedLoA, MergedLoB, MergedHiA, MergedHiB;
	MergeFwdBwd(LA, LB,
				LoA, LoB, PathFwd,
				LoA-1, LoB-1, PathBwd,
				MergedLoA, MergedLoB, MergedHiA, MergedHiB, m_XDropPath);

	return TotalScore;
	}
