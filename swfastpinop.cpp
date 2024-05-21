#include "myutils.h"
#include "mx.h"
#include "tracebit.h"
#include "xdpmem.h"

uint SWFastPinop(XDPMem &Mem,
  const int8_t * const *AP, uint LA,
  const int8_t *B, uint LB,
  int8_t Open, int8_t Ext)
	{
	int32_t *Mrow_ = myalloc(int32_t, LB+2);
	int32_t *Drow_ = myalloc(int32_t, LB+2);
	int32_t *Mrow = Mrow_ + 1;
	int32_t *Drow = Drow_ + 1;


// Use Mrow[-1], so...
	Mrow[-1] = 0;

	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = 0;
		Drow[j] = 0;
		}
	
	int32_t BestScore = 0;

// Main loop
	int32_t M0 = int32_t(0);
	for (uint i = 0; i < LA; ++i)
		{
		const int8_t *APRow = AP[i];
		int32_t I0 = 0;
		for (uint j = 0; j < LB; ++j)
			{
			int32_t SavedM0 = M0;

		// MATCH
			{
			int32_t xM = M0;
			if (Drow[j] > xM)
				xM = Drow[j];
			if (I0 > xM)
				xM = I0;
			if (0.0f >= xM)
				xM = 0;

			byte bj = B[j];
			M0 = Mrow[j];
			xM += APRow[bj];
			if (xM > BestScore)
				BestScore = xM;

			Mrow[j] = xM;
			}
			
		// DELETE
			{
			int32_t md = SavedM0 + Open;
			Drow[j] += Ext;
			if (md >= Drow[j])
				Drow[j] = md;
			}
			
		// INSERT
			{
			int32_t mi = SavedM0 + Open;
			I0 += Ext;
			if (mi >= I0)
				I0 = mi;
			}
			}
		
		M0 = 0;
		}

	myfree(Mrow_);
	myfree(Drow_);
	return BestScore;
	}
