#include "myutils.h"
#include "fragaligner.h"

static float *WeightLookup = 0;

static float GetWeight_NoLookup(float y)
	{
	float w = expf(-s_EXPFACTOR*y*y);
	return w;
	}

static float GetWeight_Lookup(float y)
	{
	int iy = int(y+0.5f);
	if (iy < 0)
		iy = 0;
	if (iy >= TBLSZ)
		iy = TBLSZ-1;
	float w2 = WeightLookup[iy];
	return w2;
	}

static bool InitWeightLookup()
	{
	WeightLookup = myalloc(float, TBLSZ);
	for (int i = 0; i < TBLSZ; ++i)
		{
		float y = float(i);
		float w = GetWeight_NoLookup(y);
		WeightLookup[i] = w;
		}
	return true;
	};
static bool InitWeightLookupDone = InitWeightLookup();

static float DALI_dpscorefun(float a, float b)
	{
	float Score = 0;
	float diff = fabs(a - b);
	float mean = (a + b)/2;
	float ratio = diff/mean;
	float w = GetWeight_Lookup(mean);
	if (mean > 100)
		Score = 0;
	else
		{
		if (mean > 0)
			Score = w*(s_DALI_d0 - ratio);
		else
			Score = w*s_DALI_d0;
		}
	return Score;
	}

float FragAligner::CalcDALIScore()
	{
	const uint LQ = m_ChainDataQ->m_L;
	const uint LF = m_ChainDataF->m_L;
	const float *DistMxVecPtrQ = m_ChainDataQ->GetDistMxPtr();
	const float *DistMxVecPtrF = m_ChainDataF->GetDistMxPtr();
	float Score = 0;
	for (uint PosF_i = 0; PosF_i < LF; ++PosF_i)
		{
		uint PosQ_i = m_LoQ + PosF_i;
		for (uint PosF_j = 0; PosF_j < LF; ++PosF_j)
			{
			uint PosQ_j = m_LoQ + PosF_j;
			float dij_Q = DistMxVecPtrQ[PosQ_i*LQ + PosQ_j];
			float dij_F = DistMxVecPtrF[PosF_i*LF + PosF_j];
			Score += DALI_dpscorefun(dij_Q, dij_F);
			}
		}
	Score += LF*s_DALI_Theta;
	return Score;
	}

void FragAligner::Align(const ChainData &Q, const ChainData &F, uint LoQ)
	{
	m_LoQ = LoQ;
	const uint LQ = m_ChainDataQ->m_L;
	const uint LF = m_ChainDataF->m_L;
	asserta(m_LoQ + LF <= LQ);
	m_DALIScore = CalcDALIScore();
	}

void cmd_test()
	{
	asserta(optset_input2);
	vector<PDBChain *> Chains1;
	vector<PDBChain *> Chains2;
	ReadChains(g_Arg1, Chains1);
	ReadChains(opt(input2), Chains2);
	
	ChainData CD_i;
	ChainData CD_j;

	FragAligner FA;
	FA.m_ChainDataQ = &CD_i;
	FA.m_ChainDataF = &CD_j;

	const uint n1 = SIZE(Chains1);
	const uint n2 = SIZE(Chains2);
	const uint FL = 32;
	for (uint i = 0; i < n1; ++i)
		{
		PDBChain &Chain_i = *Chains1[i];
		CD_i.SetChain(Chain_i);
		uint L_i = Chain_i.GetSeqLength();

		for (uint j = 0; j < n2; ++j)
			{
			PDBChain &Chain_j = *Chains2[j];
			uint L_j = Chain_j.GetSeqLength();
			CD_j.SetChain(Chain_j);
			for (uint Lo_i = 0; Lo_i < L_i; ++Lo_i)
				{
				uint Hi_i = Lo_i + L_j + 4;
				if (Hi_i >= L_i)
					break;

				FA.Align(CD_i, CD_j, Lo_i);
				Log("[%4u]  %8.3g\n", Lo_i, FA.m_DALIScore);
				}
			}
		}
	}
