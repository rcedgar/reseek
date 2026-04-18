#include "myutils.h"
#include "mx.h"
#include "tracebit.h"
#include "xdpmem.h"
#include "swtrace.h"

static inline void reverse_path_buffer(
	char *path_buffer, uint ncol)
	{
	for (uint i = 0; i < ncol/2; ++i)
		swap(path_buffer[i], path_buffer[ncol-i-1]);
	path_buffer[ncol] = 0;
	}

void TraceBackBitSW(XDPMem &Mem,
  uint LA, uint LB, uint Besti, uint Bestj,
  uint &Leni, uint &Lenj, string &Path)
	{
	Path.clear();
	byte **TB = Mem.GetTBBit();

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
				return;
				}
			t = TB[i-1][j-1];
			if (t & TRACEBITS_DM)
				State = 'D';
			else if (t & TRACEBITS_IM)
				State = 'I';
			else if (t & TRACEBITS_SM)
				{
				Leni = Besti - i + 1;
				Lenj = Bestj - j + 1;
				reverse(Path.begin(), Path.end());
				return;
				}
			else
				State = 'M';
			--i;
			--j;
			break;

		case 'D':
			asserta(i > 0);
			t = TB[i-1][j];
			if (t & TRACEBITS_MD)
				State = 'M';
			else
				State = 'D';
			--i;
			break;

		case 'I':
			asserta(j > 0);
			t = TB[i][j-1];
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
	}

float SWFast(XDPMem &Mem, const float * const *SMxData, uint LA, uint LB,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path)
	{
	asserta(Open <= 0);
	asserta(Ext <= 0);

	Mem.Clear();
	Mem.Alloc(LA+32, LB+32);
	//const float * const *SMxData = SMx.GetData();

	Leni = 0;
	Lenj = 0;

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();
	INIT_TRACE(LA, LB, TB);

#if DEBUG_UNINIT
	Mem.m_TBBit.Assign(TRACEBITS_UNINIT);
	vector<vector<bool> > TBDone(LA);
	for (uint i = 0; i < LA; ++i)
		TBDone[i].resize(LB);
#endif

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;
	TRACE_M(0, -1, MINUS_INFINITY);

	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		TRACE_M(0, j, MINUS_INFINITY);
		TRACE_D(0, j, MINUS_INFINITY);
		}
	
	float BestScore = 0.0f;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;

// Main loop
	float M0 = float (0);
	for (uint i = 0; i < LA; ++i)
		{
		const float *SMxRow = SMxData[i];
		float I0 = MINUS_INFINITY;
		byte *TBrow = TB[i];
		for (uint j = 0; j < LB; ++j)
			{
			byte TraceBits = TRACEBITS_MM;
			float SavedM0 = M0;

		// MATCH
			{
		// M0 = DPM[i][j]
		// I0 = DPI[i][j]
		// Drow[j] = DPD[i][j]
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
			TRACE_M(i, j, xM);
		// Mrow[j] = DPM[i+1][j+1])
			}
			
		// DELETE
			{
		// SavedM0 = DPM[i][j]
		// Drow[j] = DPD[i][j]
			float md = SavedM0 + Open;
			Drow[j] += Ext;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
			TRACE_D(i, j, Drow[j]);
		// Drow[j] = DPD[i+1][j]
			}
			
		// INSERT
			{
		// SavedM0 = DPM[i][j]
		// I0 = DPI[i][j]
			float mi = SavedM0 + Open;
			I0 += Ext;
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				}
			}
			
			TBrow[j] = TraceBits;
#if DEBUG_UNINIT
			TBDone[i][j] = true;
#endif
			}
		
		M0 = MINUS_INFINITY;
		}
	DONE_TRACE(BestScore, Besti, Bestj, TB);
	if (BestScore == 0.0f)
		return 0.0f;

#if DEBUG_UNINIT
	{
	for (uint i = 0; i < LA; ++i)
		{
		for (uint j = 0; j < LB; ++j)
			{
			if (TB[i][j] == TRACEBITS_UNINIT)
				Log("!TB[%u][%u]\n", i, j);
			if (!TBDone[i][j])
				Log("!TBDone[%u][%u]\n", i, j);
			}
		}
	}
#endif

	TraceBackBitSW(Mem, LA, LB, Besti+1, Bestj+1,
	  Leni, Lenj, Path);
	asserta(Besti+1 >= Leni);
	asserta(Bestj+1 >= Lenj);

	Loi = Besti + 1 - Leni;
	Loj = Bestj + 1 - Lenj;

	return BestScore;
	}

float SWFast_SubstMx(XDPMem &Mem,
	const byte *A, uint LA,
	const byte *B, uint LB,
	const vector<vector<float> > &SubstMx,
	float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
	string &Path)
	{
	asserta(Open <= 0);
	asserta(Ext <= 0);

	Mem.Clear();
	Mem.Alloc(LA+32, LB+32);
	//const float * const *SMxData = SMx.GetData();

	Leni = 0;
	Lenj = 0;

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();
	INIT_TRACE(LA, LB, TB);

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;
	TRACE_M(0, -1, MINUS_INFINITY);

	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		TRACE_M(0, j, MINUS_INFINITY);
		TRACE_D(0, j, MINUS_INFINITY);
		}
	
	float BestScore = 0.0f;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;

// Main loop
	float M0 = float (0);
	for (uint i = 0; i < LA; ++i)
		{
		byte LetterA = A[i];
		const float *SMxRow = SubstMx[LetterA].data();
		float I0 = MINUS_INFINITY;
		byte *TBrow = TB[i];
		for (uint j = 0; j < LB; ++j)
			{
			byte TraceBits = 0;
			float SavedM0 = M0;

		// MATCH
			{
		// M0 = DPM[i][j]
		// I0 = DPI[i][j]
		// Drow[j] = DPD[i][j]
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

			byte LetterB = B[j];
			M0 = Mrow[j];
			xM += SMxRow[LetterB];
			if (xM > BestScore)
				{
				BestScore = xM;
				Besti = i;
				Bestj = j;
				}

			Mrow[j] = xM;
			TRACE_M(i, j, xM);
		// Mrow[j] = DPM[i+1][j+1])
			}
			
		// DELETE
			{
		// SavedM0 = DPM[i][j]
		// Drow[j] = DPD[i][j]
			float md = SavedM0 + Open;
			Drow[j] += Ext;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
			TRACE_D(i, j, Drow[j]);
		// Drow[j] = DPD[i+1][j]
			}
			
		// INSERT
			{
		// SavedM0 = DPM[i][j]
		// I0 = DPI[i][j]
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

	DONE_TRACE(BestScore, Besti, Bestj, TB);
	if (BestScore == 0.0f)
		return 0.0f;

	TraceBackBitSW(Mem, LA, LB, Besti+1, Bestj+1,
	  Leni, Lenj, Path);
	asserta(Besti+1 >= Leni);
	asserta(Bestj+1 >= Lenj);

	Loi = Besti + 1 - Leni;
	Loj = Bestj + 1 - Lenj;

	return BestScore;
	}

using colscorefn = float(uint i, uint j);
float SWFast_Callback(XDPMem &Mem,
	uint LA, uint LB, colscorefn sf,
	float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
	string &Path)
	{
	asserta(Open <= 0);
	asserta(Ext <= 0);

	Mem.Clear();
	Mem.Alloc(LA+32, LB+32);

	Leni = 0;
	Lenj = 0;

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;

	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		TRACE_M(0, j, MINUS_INFINITY);
		TRACE_D(0, j, MINUS_INFINITY);
		}
	
	float BestScore = 0.0f;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;

// Main loop
	float M0 = float (0);
	for (uint i = 0; i < LA; ++i)
		{
		float I0 = MINUS_INFINITY;
		byte *TBrow = TB[i];
		for (uint j = 0; j < LB; ++j)
			{
			byte TraceBits = 0;
			float SavedM0 = M0;

		// MATCH
			{
		// M0 = DPM[i][j]
		// I0 = DPI[i][j]
		// Drow[j] = DPD[i][j]
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
			xM += sf(i, j);
			if (xM > BestScore)
				{
				BestScore = xM;
				Besti = i;
				Bestj = j;
				}

			Mrow[j] = xM;
			TRACE_M(i, j, xM);
		// Mrow[j] = DPM[i+1][j+1])
			}
			
		// DELETE
			{
		// SavedM0 = DPM[i][j]
		// Drow[j] = DPD[i][j]
			float md = SavedM0 + Open;
			Drow[j] += Ext;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
			TRACE_D(i, j, Drow[j]);
		// Drow[j] = DPD[i+1][j]
			}
			
		// INSERT
			{
		// SavedM0 = DPM[i][j]
		// I0 = DPI[i][j]
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

	TraceBackBitSW(Mem, LA, LB, Besti+1, Bestj+1,
	  Leni, Lenj, Path);
	asserta(Besti+1 >= Leni);
	asserta(Bestj+1 >= Lenj);

	Loi = Besti + 1 - Leni;
	Loj = Bestj + 1 - Lenj;

	return BestScore;
	}

void traceback_flat(const uint8_t *__restrict TB,
  uint LA, uint LB, uint Besti, uint Bestj,
  uint &Leni, uint &Lenj, char *path_buffer, uint &ncol)
	{
	Leni = 0;
	Lenj = 0;
	uint i = Besti;
	uint j = Bestj;
	ncol = 0;
	char State = 'M';
	for (;;)
		{
		path_buffer[ncol++] = State;

		byte t;
		switch (State)
			{
		case 'M':
			if (i == 0 || j == 0)
				{
				Leni = Besti - i;
				Lenj = Bestj - j;
				path_buffer[ncol] = 0;
				reverse_path_buffer(path_buffer, ncol);
				return;
				}
			//t = TB[i-1][j-1];
			t = TB[(i-1)*LB + j-1];
			if (t & TRACEBITS_DM)
				State = 'D';
			else if (t & TRACEBITS_IM)
				State = 'I';
			else if (t & TRACEBITS_SM)
				{
				Leni = Besti - i + 1;
				Lenj = Bestj - j + 1;
				path_buffer[ncol] = 0;
				reverse_path_buffer(path_buffer, ncol);
				return;
				}
			else
				State = 'M';
			--i;
			--j;
			break;

		case 'D':
			asserta(i > 0);
			//t = TB[i-1][j];
			t = TB[(i-1)*LB + j];
			if (t & TRACEBITS_MD)
				State = 'M';
			else
				State = 'D';
			--i;
			break;

		case 'I':
			asserta(j > 0);
			// t = TB[i][j-1];
			t = TB[i*LB + j-1];
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
	}

// scratch_rows length 2*LB + 2
// TB length LA*LB
float sw_flat(
	float *__restrict scratch_rows,
	uint8_t *__restrict TB,
	uint LA, uint LB, colscorefn sf,
	float Open, float Ext, uint &Loi, uint &Loj,
	char *path_buffer, uint &ncol)
	{
	asserta(Open <= 0);
	asserta(Ext <= 0);

	memset(TB, 0, LA*LB);//@@TODO

	float *Mrow = scratch_rows + 1;
	float *Drow = scratch_rows + LB + 2;

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;

	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		TRACE_M(0, j, MINUS_INFINITY);
		TRACE_D(0, j, MINUS_INFINITY);
		}
	
	float BestScore = 0.0f;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;

// Main loop
	float M0 = float(0);
	for (uint i = 0; i < LA; ++i)
		{
		float I0 = MINUS_INFINITY;
		uint8_t *TBrow = TB + i*LB;
		for (uint j = 0; j < LB; ++j)
			{
			byte TraceBits = 0;
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
			xM += sf(i, j);
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

	uint Leni = UINT_MAX;
	uint Lenj = UINT_MAX;
	traceback_flat(TB, LA, LB, Besti+1, Bestj+1,
	  Leni, Lenj, path_buffer, ncol);
	asserta(Besti+1 >= Leni);
	asserta(Bestj+1 >= Lenj);

	Loi = Besti + 1 - Leni;
	Loj = Bestj + 1 - Lenj;

	return BestScore;
	}

// scratch_rows length 2*LT + 2
// TB length LQ*LT
// scratch_ppsms length nfeat
float sw_flat_pssm(
	float *__restrict scratch_rows,
	uint8_t *__restrict TB,
	const float ** __restrict scratch_ppsms,
	const uint8_t *__restrict profQ, uint LQ,
	const float *__restrict pssmT, uint LT,
	const uint32_t * __restrict feature_block_offsets,
	uint nfeat,
	float Open, float Ext,
	uint &loQ, uint &loT,
	char *path_buffer, uint &ncol)
	{
	asserta(Open <= 0);
	asserta(Ext <= 0);

	float *Mrow = scratch_rows + 1;
	float *Drow = scratch_rows + LT + 2;

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;

	for (uint j = 0; j <= LT; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		TRACE_M(0, j, MINUS_INFINITY);
		TRACE_D(0, j, MINUS_INFINITY);
		}

	float BestScore = 0.0f;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;

// Main loop
	float M0 = float(0);
	for (uint i = 0; i < LQ; ++i)
		{
		// Select one PSSM row per feature for this i.
		for (uint fi = 0; fi < nfeat; ++fi)
			{
			const uint8_t codeA = profQ[size_t(fi)*LQ + i];
			const float * __restrict pssm_fi =
				pssmT + size_t(feature_block_offsets[fi])*LT;
			scratch_ppsms[fi] = pssm_fi + size_t(codeA)*LT;
			}

		float I0 = MINUS_INFINITY;
		uint8_t *TBrow = TB + i*LT;
		for (uint j = 0; j < LT; ++j)
			{
			byte TraceBits = 0;
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

			float Score = 0.0f;
			for (uint fi = 0; fi < nfeat; ++fi)
				Score += scratch_ppsms[fi][j];

			xM += Score;

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

	uint Leni = UINT_MAX;
	uint Lenj = UINT_MAX;
	traceback_flat(TB, LQ, LT, Besti+1, Bestj+1, Leni, Lenj,
		path_buffer, ncol);

	assert(Leni <= Besti+1);
	assert(Lenj <= Bestj+1);

	loQ = Besti + 1 - Leni;
	loT = Bestj + 1 - Lenj;

	assert(loQ < LQ);
	assert(loT < LT);

	return BestScore;
	}
