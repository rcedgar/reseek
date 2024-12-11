#include "myutils.h"
#include "mx.h"
#include "xdpmem.h"
#if 0
float GaplessSWFast(XDPMem &Mem, const Mx<float> &SMx, uint LA, uint LB)
	{
#if TRACE && !DOTONLY
	SMx.LogMe();
#endif

	Mem.Alloc(LA+32, LB+32);
	const float * const *SMxData = SMx.GetData();

	float *Mrow = Mem.GetDPRow1();

	Mrow[-1] = MINUS_INFINITY;
	for (uint j = 0; j <= LB; ++j)
		Mrow[j] = MINUS_INFINITY;
	
	float BestScore = 0.0f;
	float M0 = float(0);
	for (uint i = 0; i < LA; ++i)
		{
		const float *SMxRow = SMxData[i];
		float I0 = MINUS_INFINITY;
		for (uint j = 0; j < LB; ++j)
			{
			float xM = M0;
			if (0.0f >= xM)
				xM = 0.0f;

			M0 = Mrow[j];
			xM += SMxRow[j];
			if (xM > BestScore)
				BestScore = xM;

			Mrow[j] = xM;
			}
		M0 = MINUS_INFINITY;
		}
	return BestScore;
	}
#endif