#if 0
#include "myutils.h"
#include "tracebit.h"

static const float MINUS_INFINITY = std::numeric_limits<float>::lowest();

// scratch_float size 2*LB + 2
// scratch_byte size LA*LB
float sw_flat(
	const float * __restrict const smx,
	float * __restrict scratch_float,
	uint8_t * __restrict scratch_byte,
	uint LA,
	uint LB,
	float Open,
	float Ext,
	uint &Loi,
	uint &Loj,
	uint &Leni,
	uint &Lenj,
	string &Path)
	{
	asserta(Open <= 0);
	asserta(Ext <= 0);

	Leni = 0;
	Lenj = 0;

	float *Mrow = scratch_float + 1;
	float *Drow = scratch_float + LB + 2;

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;

	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}
	
	float BestScore = 0.0f;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;
	float M0 = float (0);
	for (uint i = 0; i < LA; ++i)
		{
		const float * __restrict SMxRow = smx + LB*i;
		float I0 = MINUS_INFINITY;
		byte *__restrict TBrow = scratch_byte + LB*i;
		for (uint j = 0; j < LB; ++j)
			{
			byte TraceBits = TRACEBITS_MM;
			float SavedM0 = M0;

		// MATCH
			{
			float xM = M0;
			if (Drow[j] > xM)
				{
				xM = Drow[j];
				TraceBits = TRACEBITS_DM;
				}
			if (I0 > xM)
				{
				xM = I0;
				TraceBits = TRACEBITS_IM;
				}
			if (0.0f >= xM)
				{
				xM = 0.0f;
				TraceBits = TRACEBITS_SM;
				}

			M0 = Mrow[j];
			xM += SMxRow[j];
			if (xM > BestScore)
				{
				BestScore = xM;
				Besti = i;
				Bestj = j;
				}

			Mrow[j] = xM;
			}
			
		// DELETE
			{
			float md = SavedM0 + Open;
			Drow[j] += Ext;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
			}
			
		// INSERT
			{
			float mi = SavedM0 + Open;
			I0 += Ext;
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				}
			}
			
			TBrow[j] = TraceBits;
			}
		
		M0 = MINUS_INFINITY;
		}
	if (BestScore == 0.0f)
		return 0.0f;

	{ // Traceback
	Path.clear();
	Path.reserve(2*max(LA,LB));
	byte * __restrict TB = scratch_byte;

	uint i = Besti;
	uint j = Bestj;
	char State = 'M';
	for (;;)
		{
		Path += State;

		byte t;
		switch (State)
			{
		case 'M':
			if (i == 0 || j == 0)
				{
				Leni = Besti - i;
				Lenj = Bestj - j;
				reverse(Path.begin(), Path.end());
				goto Done;
				}
			t = TB[(i-1)*LB + j - 1];
			if (t & TRACEBITS_DM)
				State = 'D';
			else if (t & TRACEBITS_IM)
				State = 'I';
			else if (t & TRACEBITS_SM)
				{
				Leni = Besti - i + 1;
				Lenj = Bestj - j + 1;
				reverse(Path.begin(), Path.end());
				goto Done;
				}
			else
				State = 'M';
			--i;
			--j;
			break;

		case 'D':
			asserta(i > 0);
			t = TB[(i-1)*LB + j];
			if (t & TRACEBITS_MD)
				State = 'M';
			else
				State = 'D';
			--i;
			break;

		case 'I':
			asserta(j > 0);
			t = TB[i*LB + j - 1];
			if (t & TRACEBITS_MI)
				State = 'M';
			else
				State = 'I';
			--j;
			break;

		default:
			Die("TraceBackBitSW, invalid state %c", State);
			}
		}
	Done:;
	}
	asserta(Besti+1 >= Leni);
	asserta(Bestj+1 >= Lenj);

	Loi = Besti + 1 - Leni;
	Loj = Bestj + 1 - Lenj;

	return BestScore;
	}
#endif // 0