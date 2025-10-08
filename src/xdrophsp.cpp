#include "myutils.h"
#include "dssaligner.h"

#if TRACE_XDROP
void LogAln(const char *A, const char *B, const char *Path, unsigned ColCount);
#endif

float DSSAligner::SubstScore(uint PosA, uint PosB)
	{
	const vector<vector<byte> > &ProfileA = *m_ProfileA;
	const vector<vector<byte> > &ProfileB = *m_ProfileB;
	const uint FeatureCount = DSSParams::GetFeatureCount();
	assert(SIZE(ProfileA) == FeatureCount);
	assert(SIZE(ProfileB) == FeatureCount);
	float Total = 0;
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = DSSParams::m_Features[FeatureIdx];
		uint AlphaSize = g_AlphaSizes2[F];
		float **ScoreMx = DSSParams::m_ScoreMxs[F];
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

float DSSAligner::XDropHSP(uint Loi_in, uint Loj_in, uint Len,
						   uint &Loi_out, uint &Loj_out,
						   uint &Hii_out, uint &Hij_out)
	{
#if TRACE_XDROP
	{
	Log("\nXDropHSP(Loi=%u, Loj=%u, Len=%u\n", Loi_in, Loj_in, Len);
	string PathHSP;
	for (uint i = 0; i < Len; ++i)
		PathHSP += 'M';
	LogAln(m_ChainA->m_Seq.c_str() + Loi_in, m_ChainB->m_Seq.c_str() + Loj_in,
		   PathHSP.c_str(), Len);
	}
#endif
	Loi_out = UINT_MAX;
	Loj_out = UINT_MAX;
	Hii_out = UINT_MAX;
	Hij_out = UINT_MAX;
	const float Open = DSSParams::m_GapOpen;
	const float Ext = DSSParams::m_GapExt;
	const float X = float(DSSParams::m_MKF_X2);
	const uint LA = m_ChainA->GetSeqLength();
	const uint LB = m_ChainB->GetSeqLength();

	uint LoA = Loi_in + Len/2;
	uint LoB = Loj_in + Len/2;

// Find highest-scoring 8mer
	const uint K = 8;
	asserta(Len > K);

	float *v = myalloc(float, Len);
	for (uint Col = 0; Col < Len; ++Col)
		v[Col] = StaticSubstScore((void *) this, Loi_in + Col, Loj_in + Col);

	float BestMerScore = 0;
	for (uint MerStart = 0; MerStart + K <= Len; ++MerStart)
		{
		float MerScore = 0;
		for (uint k = 0; k < K; ++k)
			MerScore += v[MerStart+k];
		if (MerScore > BestMerScore)
			{
			BestMerScore = MerScore;
			LoA = Loi_in + MerStart;
			LoB = Loj_in + MerStart;
			}
		}
	myfree(v);

	//asserta(LoA > 0 && LoB > 0);
	uint MinLo = min(LoA, LoB);
	if (MinLo < K/2)
		{
		LoA += K/2;
		LoB += K/2;
		}

	string FwdPath, BwdPath;
	uint FwdSegLoA, FwdSegLoB;
	float ScoreFwd = 
		XDropFwd(m_Mem, X, Open, Ext, StaticSubstScore, (void *) this, 
				 LoA, LA, LoB, LB, &FwdSegLoA, &FwdSegLoB, FwdPath);

	uint BwdSegLoA, BwdSegLoB;
	float ScoreBwd = 
		XDropBwd(m_Mem, X, Open, Ext, StaticSubstScore, (void *) this, 
				 LoA-1, LA, LoB-1, LB, &BwdSegLoA, &BwdSegLoB, BwdPath);
	float TotalScore = ScoreFwd + ScoreBwd;
	if (TotalScore < 10)
		{
		m_XDropPath.clear();
		return 0;
		}

	MergeFwdBwd(LA, LB,
				LoA, LoB, FwdPath,
				LoA-1, LoB-1, BwdPath,
				Loi_out, Loj_out, Hii_out, Hij_out, m_XDropPath);

#if TRACE_XDROP
	{
	uint BwdM, BwdD, BwdI;
	void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);
	GetPathCounts(BwdPath, BwdM, BwdD, BwdI);

	uint BwdColsA = BwdM + BwdD;
	uint BwdColsB = BwdM + BwdI;
	asserta(LoA >= BwdColsA);
	asserta(LoB >= BwdColsB);
	uint BwdLoA = LoA - BwdColsA;
	uint BwdLoB = LoB - BwdColsB;

	Log("\nFwd:\n");
	LogAln(m_ChainA->m_Seq.c_str() + LoA, m_ChainB->m_Seq.c_str() + LoB,
		   FwdPath.c_str(), SIZE(FwdPath));

	Log("\nBwd:\n");
	LogAln(m_ChainA->m_Seq.c_str() + BwdLoA, m_ChainB->m_Seq.c_str() + BwdLoB,
		   BwdPath.c_str(), SIZE(BwdPath));

	Log("\nMerged:\n");
	LogAln(m_ChainA->m_Seq.c_str() + Loi_out, m_ChainB->m_Seq.c_str() + Loj_out,
		   m_XDropPath.c_str(), SIZE(m_XDropPath));
	}
#endif

	return TotalScore;
	}
