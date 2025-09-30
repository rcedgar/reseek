#include "myutils.h"
#include "mx.h"
#include "alpha.h"
#include "xdpmem.h"

float SWFastGaplessProfb(float *DProw_, const float * const *ProfA, uint LA,
  const byte *B, uint LB)
	{
	float *MrowF = DProw_ + 1;
	float *MrowR = DProw_ + LB + 2;

// Use Mrow[-1], so...
	for (uint j = 0; j <= 2*LB + 3; ++j)
		DProw_[j] = MINUS_INFINITY;
	
	float BestScoreF = 0.0f;
	float BestScoreR = 0.0f;

// Main loop
	float M0F = float (0);
	float M0R = float (0);
	for (uint i = 0; i < LA; ++i)
		{
		const float *SMxRowF = ProfA[i];
		const float *SMxRowR = ProfA[LA-i-1];
		for (uint j = 0; j < LB; ++j)
			{
			float SavedM0F = M0F;
			float SavedM0R = M0R;

		// MATCH
			{
			float xMF = M0F;
			float xMR = M0R;
			if (xMF < 0.0f)
				xMF = 0.0f;
			if (xMR < 0.0f)
				xMR = 0.0f;

			M0F = MrowF[j];
			M0R = MrowR[j];
			byte bj = B[j];
			assert(bj < 36);
			xMF += SMxRowF[bj];
			xMR += SMxRowR[bj];
			if (xMF > BestScoreF)
				BestScoreF = xMF;
			if (xMR > BestScoreR)
				BestScoreR = xMR;

			MrowF[j] = xMF;
			MrowR[j] = xMR;
			}
			}
		
		M0F = 0;
		M0R = 0;
		}
	float Score = BestScoreF - BestScoreR;
	return Score;
	}
