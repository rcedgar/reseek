#pragma once

#include "pdbchain.h"
#include "chaindata.h"
#include "maxbuff.h"

static const uint DEFAULT_MAX_LF = 64;
static const float s_DALI_D = 20.0f;
static const float s_DALI_d0 = 0.2f;
static const float s_DALI_Theta = 0.2f;
static const float s_EXPFACTOR = 1.0f/(s_DALI_D*s_DALI_D);
static const int TBLSZ = 100;

class FragAligner
	{
public:
	ChainData *m_ChainDataQ = 0;
	ChainData *m_ChainDataF = 0;
	const float *m_SubstMxPtr = 0;
	uint m_LoQ = UINT_MAX;
	float m_DALIScore = FLT_MAX;

public:
	FragAligner() {}
	~FragAligner() {}

	void Align(const ChainData &A, const ChainData &F, uint LoQ);

private:
	float CalcDALIScore();
	};
