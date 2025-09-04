#include "myutils.h"
#include "mx.h"
#include "xdpmem.h"
#include "tracebit.h"
#include "pathinfo.h"
#include "swtrace.h"

void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);

static void TraceBack(XDPMem &Mem, uint Besti, uint Bestj, string &Path)
	{
	Path.clear();
	byte **TB = Mem.GetTBBit();
	uint i = Besti;
	uint j = Bestj;
	char State = 'M';
	for (;;)
		{
		Path += State;
#if TRACE && !DOTONLY
		Log("Traceback %u, %u, %c  \"%s\"\n", i, j, State, Path.c_str());
#endif
		//if (i == 1 && j == 1)
		if (i == 1 || j == 1)
			break;

		char NextState = '?';
		switch (State)
			{
		case 'M':
			NextState = GetTBBitM(TB, i, j);
#if TRACE && !DOTONLY
			Log(" GetTBBitM(%u, %u) = %c\n", i, j, NextState);
#endif
			asserta(i > 0 && j > 0);
			--i;
			--j;
			break;

		case 'D':
			NextState = GetTBBitD(TB, i, j);
#if TRACE && !DOTONLY
			Log(" GetTBBitD(%u, %u) = %c\n", i, j, NextState);
#endif
			asserta(i > 0);
			--i;
			break;

		case 'I':
			NextState = GetTBBitI(TB, i, j);
#if TRACE && !DOTONLY
			Log(" GetTBBitI(%u, %u) = %c\n", i, j, NextState);
#endif
			asserta(j > 0);
			--j;
			break;

		default:
			Die("TraceBackBit, invalid state %c", State);
			}
		State = NextState;
		}
	std::reverse(Path.begin(), Path.end());
#if TRACE && !DOTONLY
	Log("Traceback = %s\n", Path.c_str());
#endif
	}

// After traceback start of path is *ptrSegLoA,*ptrSegLoB
// may be different from LoA,LoB
float XDropFwd(XDPMem &Mem,
  float X, float Open, float Ext, 
  fn_SubstScore SubFn, void *UserData,
  uint LoA, uint aLA, uint LoB, uint aLB,
  uint *ptrSegLoA, uint *ptrSegLoB,
  string &Path)
	{
	*ptrSegLoA = UINT_MAX;
	*ptrSegLoB = UINT_MAX;
	asserta(LoA < aLA);
	asserta(LoB < aLB);
	uint LA = aLA - LoA;
	uint LB = aLB - LoB;

	Path.clear();

	if (LA == 1 || LB == 1)
		{
		float Score = SubFn(UserData, LoA, LoB);
		if (Score > 0)
			Path.push_back('M');
		return Score;
		}

	if (Open > 0.0f || Ext > 0.0f)
		Warning("XDropFwdFast(): non-negative open %.1f, ext %.1f", Open, Ext);
	const float AbsOpen = -Open;
	const float AbsExt = -Ext;

	Mem.Alloc(LA+1, LB+1);

	byte **TB = Mem.GetTBBit();
	INIT_TRACE(LA, LB, TB);

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();

	Mrow[-1] = MINUS_INFINITY;
	TRACE_M(0, -1, MINUS_INFINITY);

	Drow[0] = MINUS_INFINITY;
	Drow[1] = MINUS_INFINITY;
	TRACE_D(0, 0, MINUS_INFINITY);
	TRACE_D(0, 1, MINUS_INFINITY);

// Main loop
	float BestScore = 0;
	TRACE_M(1, 1, 0);
	uint Besti = 0;
	uint Bestj = 0;

	uint prev_jlo = 0;
	uint prev_jhi = 0;
	uint jlo = 1;
	uint jhi = 1;

// Inner loop does this:
//	Mrow[j] = DPM[i][j+1] -> DPM[i+1][j+1]
//	Drow[j] = DPD[i][j]   -> DPD[i+1][j]

	float M0 = BestScore;
	for (uint i = 1; i <= LA; ++i)
		{
#if TRACE && !DOTONLY
		Log("XDrop i=%u j=%u .. %u\n", i, jlo, jhi);
#endif
		if (jlo == prev_jlo)
			{
			asserta(jlo>0);
			Mrow[jlo-1] = MINUS_INFINITY;
			Drow[jlo] = MINUS_INFINITY;
			TRACE_M(i, jlo-1, MINUS_INFINITY);
			TRACE_D(i, jlo, MINUS_INFINITY);
			}

		uint endj = min(prev_jhi+1,LB);
		for (uint j = endj+1; j <= min(jhi+1, LB); ++j)
			{
			Mrow[j-1] = MINUS_INFINITY;
			Drow[j] = MINUS_INFINITY;
			TRACE_M(i, j-1, MINUS_INFINITY);
			TRACE_D(i, j, MINUS_INFINITY);
			}

		uint next_jlo = UINT_MAX;
		uint next_jhi = UINT_MAX;

		float I0 = MINUS_INFINITY;

		byte *TBrow = TB[i];
		asserta(jlo>0);
		asserta(jlo<=jhi);
		float SavedM0 = UNINIT;

		for (uint j = jlo; j <= jhi; ++j)
			{
			byte TraceBits = 0;

			SavedM0 = M0; // SavedM0 = M0 = DPM[i][j]

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
			M0 = Mrow[j]; // M0 = DPM[i][j+1]

			//float s = xM + MxRow[b];
			float s = SubFn(UserData, LoA + i-1, LoB + j-1);
			TRACE_Sub(i-1, j-1, s);
			s += xM;
			Mrow[j] = s;	// Mrow[j] = DPM[i+1][j+1]
			TRACE_M(i, j, s);

			float h = s - BestScore + X;
		// Match-Match
			if (h > 0)
				{
				next_jlo = min(next_jlo, j+1);
				next_jhi = j+1;
				}

		// Match-Delete
			if (h > AbsOpen)
				next_jlo = min(next_jlo, j);

		// Match-Insert potentially extends current row
			if (h > AbsExt && j == jhi && jhi + 1 < LB)
				{
				++jhi;
				uint new_endj = min(jhi+1, LB);
				new_endj = max(new_endj, endj);
				for (uint j2 = endj+1; j2 <= new_endj; ++j2)
					{
				// Nasty special case for j=j2-1, Mrow[j] has already
				// been updated for current i.
					if (j2-1 > j)
						{
						Mrow[j2-1] = MINUS_INFINITY;
						TRACE_M(i, j2-1, MINUS_INFINITY);
						}

					Drow[j2] = MINUS_INFINITY;
					TRACE_M(i, j2, MINUS_INFINITY);
					}
				endj = new_endj;
				}
			if (s >= BestScore)
				{
				BestScore = s;
				Besti = i;
				Bestj = j;
				}
			}

		// DELETE
			if (j != jlo)
				{
		// SavedM0 = DPM[i][j]
		// Drow[j] = DPD[i][j]
			float md = SavedM0 + Open;
			Drow[j] += Ext;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				TRACE_D(i, j, md);
				}
		// Drow[j] = DPD[i+1][j]
			float h = Drow[j] - BestScore + X;

		// Delete-Match
			if (h > 0)
				{
				next_jlo = min(next_jlo, j-1);
				next_jhi = max(next_jhi, j-1);
				}
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
		// I0 = DPI[i][j+1]

			float h = I0 - BestScore + X;
		// Insert-Match
			if (h > 0)
				{
				next_jlo = min(next_jlo, j+1);
				next_jhi = max(next_jhi, j+1);
				}

		// Insert-Insert potentially extends current row
			if (h > AbsExt && j == jhi && jhi + 1 < LB)
				{
				++jhi;
				uint new_endj = min(jhi+1, LB);
				new_endj = max(new_endj, endj);
				for (uint j2 = endj+1; j2 <= new_endj; ++j2)
					{
					Mrow[j2-1] = MINUS_INFINITY;
					Drow[j2] = MINUS_INFINITY;
					TRACE_M(i, j2-1, MINUS_INFINITY);
					TRACE_D(i, j, MINUS_INFINITY);
					}
				endj = new_endj;
				}
			}
		
			TBrow[j] = TraceBits;
			}

	// Special case for end of Drow[]
		if (jhi < LB)
			{
			const uint jhi1 = jhi+1;
			TBrow[jhi1] = 0;
			float md = M0 + Open;
			Drow[jhi1] += Ext;
			if (md >= Drow[jhi1])
				{
				Drow[jhi1] = md;
				TBrow[jhi1] = TRACEBITS_MD;
				TRACE_D(i, jhi1, md);
				}
			}

		if (next_jlo == UINT_MAX)
			break;

		prev_jlo = jlo;
		prev_jhi = jhi;
		jlo = next_jlo;
		jhi = next_jhi;
		if (jlo > LB)
			jlo = LB;
		if (jhi > LB)
			jhi = LB;
		asserta(jlo <= jhi);
		asserta(jlo >= prev_jlo);

		if (jlo == prev_jlo)
			{
			M0 = MINUS_INFINITY;
			Drow[jlo] = MINUS_INFINITY;
			TRACE_D(i, jlo, MINUS_INFINITY);
			}
		else
			{
			assert(jlo > prev_jlo);
			M0 = Mrow[jlo-1];
			}
		}

	DONE_TRACE(BestScore, Besti, Bestj, TB);
	if (BestScore <= 0.0f)
		return 0.0f;

	TraceBack(Mem, Besti, Bestj, Path);
	uint nM, nD, nI;
	GetPathCounts(Path, nM, nD, nI);
	uint Loi = LoA + Besti - nM - nD;
	uint Loj = LoB + Bestj - nM - nI;
	*ptrSegLoA = Loi;
	*ptrSegLoB = Loj;

#if DEBUG
	{
	const uint ColCount = SIZE(Path);
	uint PosA = Loi;
	uint PosB = Loj;
	float Score2 = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M': 
			Score2 += SubFn(UserData, PosA, PosB);
			++PosA;
			++PosB;
			break;

		case 'D':
			++PosA;
			break;

		case 'I':
			++PosB;
			break;
			}
		}
	asserta(Score2 + 0.1 >= BestScore);
	}
#endif
	return BestScore;
	}
