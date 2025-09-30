#if 0 // @@DELETE
#include "myutils.h"
#include "mx.h"
#include "tracebit.h"
#include "xdpmem.h"

uint SWFastPinopGapless(const int8_t * const *AP, uint LA,
  const int8_t *B, uint LB)
	{
	int32_t *Mrow_ = myalloc(int32_t, LB+2);
	int32_t *Mrow = Mrow_ + 1;

// Use Mrow[-1], so...
	Mrow[-1] = 0;
	for (uint j = 0; j <= LB; ++j)
		Mrow[j] = 0;
	
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

			int32_t xM = M0;
			if (xM < 0.0f)
				xM = 0;

			byte bj = B[j];
			M0 = Mrow[j];
			xM += APRow[bj];
			if (xM > BestScore)
				BestScore = xM;

			Mrow[j] = xM;
			}
		
		M0 = 0;
		}

	myfree(Mrow_);
	return BestScore;
	}
#endif