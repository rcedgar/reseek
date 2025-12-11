#include "myutils.h"
#include "fragaligner.h"

int FragAligner::m_DistScale = 9;

float FragAligner::Align(const uint *DistsPtrA, const uint *DistsPtrB) const
	{
	float Score = 0;
	for (uint i = 0; i < m_IdxCount; ++i)
		{
		uint dij_A = DistsPtrA[i];
		uint dij_B = DistsPtrB[i];
		assert(dij_A < m_DistN);
		assert(dij_B < m_DistN);
		Score += m_DistPairScoresPtr[m_DistN*dij_A + dij_B];
		}
	return Score;
	}

void FragAligner::GetIdx_ijs(uint Length, uint BandWidth,
	vector<uint> &Idx_is, vector<uint> &Idx_js)
	{
	Idx_is.clear();
	Idx_js.clear();
	for (uint i = 0; i < Length; ++i)
		{
		for (uint j = i + 1 + BandWidth; j < Length; ++j)
			{
			Idx_is.push_back(i);
			Idx_js.push_back(j);
			}
		}
	}

float *FragAligner::MakeDistPairScoresPtr(uint DistN)
	{
	float *DistPairScoresPtr = myalloc(float, DistN*DistN);
	for (uint Dist1 = 0; Dist1 < DistN; ++Dist1)
		for (uint Dist2 = 0; Dist2 < DistN; ++Dist2)
			DistPairScoresPtr[Dist1*DistN + Dist2] =
				GetDistPairScore(Dist1, Dist2);
	return DistPairScoresPtr;
	}

void FragAligner::Init(uint Length, uint BandWidth, uint DistN,
	const float *DistPairScoresPtr)
	{
	m_Length = Length;
	m_BandWidth = BandWidth;
	m_DistN = DistN;
	m_DistPairScoresPtr = DistPairScoresPtr;
	GetIdx_ijs(Length, BandWidth, m_Idx_is, m_Idx_js);
	m_IdxCount = SIZE(m_Idx_is);
	asserta(SIZE(m_Idx_js) == m_IdxCount);
	}

static float Weight_NoLookup(float y)
	{
	const float D = 20.0;
	const float x = 1/(D*D);
	float w = expf(-x*y*y);
	return w;
	}

float FragAligner::GetDistPairScore(uint Dist1, uint Dist2)
	{
	const float d0 = 0.2f;

	float a = float(Dist1);
	float b = float(Dist2);
	float Score = 0;
	float diff = fabs(a - b);
	float mean = (a + b)/2;
	float ratio = diff/mean;
	float w = Weight_NoLookup(mean);
	if (mean > 100)
		Score = 0;
	else if (mean > 0)
		Score = w*(d0 - ratio);
	else
		Score = w*d0;
	return Score;
	}

uint *FragAligner::GetDistsPtr(const PDBChain &Chain, uint Pos) const
	{
	const uint L = Chain.GetSeqLength();
	const uint n = SIZE(m_Idx_is);
	asserta(n > 0);
	asserta(SIZE(m_Idx_js) == n);
	uint *DistsPtr = myalloc(uint, n);
	for (uint k = 0; k < n; ++k)
		{
		uint i = m_Idx_is[k];
		uint j = m_Idx_js[k];
		double d = Chain.GetDist(Pos+i, Pos+j);
		uint ud = uint(d);
		if (ud >= m_DistN)
			ud = m_DistN - 1;
		DistsPtr[k] = ud;
		}
	return DistsPtr;
	}

void FragAligner::LogScoreMx() const
	{
	Log("FragAligner::LogScoreMx()\n");
	Log("Length  %u\n", m_Length);
	Log("Band    %u\n", m_BandWidth);
	Log("DistN   %u\n", m_DistN);
	Log("NrIdxs  %u\n", m_IdxCount);
	Log("\n");
	Log("Idx_ijs\n");
	for (uint i = 0; i < m_IdxCount; ++i)
		Log("  %3u  %3u\n", m_Idx_is[i], m_Idx_js[i]);

	for (uint Dist1 = 0; Dist1 < m_DistN; ++Dist1)
		{
		for (uint Dist2 = Dist1; Dist2 < m_DistN; ++Dist2)
			{
			float Score = m_DistPairScoresPtr[Dist1*m_DistN + Dist2];
			float Score2 = GetDistPairScore(Dist1, Dist2);
			Log("  %3u  %3u  %7.3f\n", Dist1, Dist2, Score);
			asserta(Score == Score2);
			}
		}
	}

#if 0
void cmd_test()
	{
	asserta(optset_input2);
	vector<PDBChain *> Chains;
	ReadChains(g_Arg1, Chains);
	const uint ChainCount = SIZE(Chains);
	ProgressLog("%u chains\n", ChainCount);

	const uint Length = 32;
	const uint BandWidth = 2;
	const uint DistN = 20;
	float *DistPairScoresPtr = FragAligner::MakeDistPairScoresPtr(DistN);

	FragAligner FA;
	FA.Init(Length, BandWidth, DistN, DistPairScoresPtr);
	//FA.LogScoreMx();

	float MinDiagScore = FLT_MAX;
	float MaxDiagScore = FLT_MAX;
	float MinOtherScore = FLT_MAX;
	float MaxOtherScore = FLT_MAX;
	for (uint i = 0; i < ChainCount; ++i)
		{
		PDBChain &Chain = *Chains[i];
		const uint L = Chain.GetSeqLength();
		for (uint Pos = 0; Pos + Length < L; ++Pos)
			{
			uint *DistsPtr = FA.GetDistsPtr(Chain, Pos);
			float Score = FA.Align(DistsPtr, DistsPtr);
			if (MinDiagScore == FLT_MAX)
				MinDiagScore = MaxDiagScore = Score;
			else
				{
				MinDiagScore = min(Score, MinDiagScore);
				MaxDiagScore = max(Score, MaxDiagScore);
				}
			}
		for (uint Pos1 = 0; Pos1 + Length < L; ++Pos1)
			{
			uint *DistsPtr1 = FA.GetDistsPtr(Chain, Pos1);
			for (uint Pos2 = Pos1 + 1; Pos2 + Length < L; ++Pos2)
				{
				uint *DistsPtr2 = FA.GetDistsPtr(Chain, Pos2);
				float Score = FA.Align(DistsPtr1, DistsPtr2);
				if (MinOtherScore == FLT_MAX)
					MinOtherScore = MaxOtherScore = Score;
				else
					{
					MinOtherScore = min(Score, MinOtherScore);
					MaxOtherScore = max(Score, MaxOtherScore);
					}
				}
			}
		}
	ProgressLog(" Diag scores %8.3g to %8.3g\n",
		MinDiagScore, MaxDiagScore);
	ProgressLog("Other scores %8.3g to %8.3g\n",
		MinOtherScore, MaxOtherScore);
	}
#endif