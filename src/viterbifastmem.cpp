#include "myutils.h"
#include "xdpmem.h"
#include "tracebit.h"
#include "pathinfo.h"

static float s_Open = -1;
static float s_Ext = -0.05f;
static float s_TermOpen = 0;
static float s_TermExt = 0;

void TraceBackBitMem(XDPMem &Mem, uint LA, uint LB, char State, string &Path);

void InitGapStr()
	{
	asserta(optset_gapstr);
	vector<string> Fields;
	Split(opt_gapstr, Fields, '_');
	asserta(SIZE(Fields) == 4);
	s_Open = StrToFloatf(Fields[0]);
	s_Ext = StrToFloatf(Fields[1]);
	s_TermOpen = StrToFloatf(Fields[2]);
	s_TermExt = StrToFloatf(Fields[3]);
	ProgressLog("open=%.3g, ext=%.3g, term_open=%.3g, term_ext=%.3g\n",
				s_Open, s_Ext, s_TermOpen, s_TermExt);

	s_Open *= -1;
	s_Ext *= -1;
	s_TermOpen *= -1;
	s_TermExt *= -1;
	}

float ViterbiFastMem(XDPMem &Mem, uint LA, uint LB,
	fn_SubstScore SubFn, void *UserData, string &Path)
	{
	if (LA*LB > 100*1000*1000)
		Die("ViterbiFastMem, seqs too long LA=%u, LB=%u", LA, LB);

	Mem.Clear();
	Mem.Alloc(LA+32, LB+32);
	
	float Open = s_TermOpen;
	float Ext = s_TermExt;

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;
	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}
	
// Main loop
	float M0 = float(0);
	for (uint i = 0; i < LA; ++i)
		{
		//const float *MxRow = Mx[a];
		float Open = s_TermOpen;
		float Ext = s_TermExt;
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
			M0 = Mrow[j];

			float Sub = SubFn(UserData, i, j);
			Mrow[j] = xM + Sub;
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
		// I0 = DPI[i][j+1]
			}
			
			Open = s_Open;
			Ext = s_Ext;
			
			TBrow[j] = TraceBits;
			}
		
	// Special case for end of Drow[]
		{
	// M0 = DPM[i][LB]
	// Drow[LB] = DPD[i][LB]
		
		TBrow[LB] = 0;
		float md = M0 + s_TermOpen;
		Drow[LB] += s_TermExt;
		if (md >= Drow[LB])
			{
			Drow[LB] = md;
			TBrow[LB] = TRACEBITS_MD;
			}
	// Drow[LB] = DPD[i+1][LB]
		}
		
		M0 = MINUS_INFINITY;

		Open = s_Open;
		Ext = s_Ext;
		}
	
// Special case for last row of DPI
	byte *TBrow = TB[LA];
	float I1 = MINUS_INFINITY;
	for (uint j = 1; j < LB; ++j)
		{
	// Mrow[j-1] = DPM[LA][j]
	// I1 = DPI[LA][j]
		
		TBrow[j] = 0;
		float mi = Mrow[int(j)-1] + s_TermOpen;
		I1 += s_TermExt;
		if (mi > I1)
			{
			I1 = mi;
			TBrow[j] = TRACEBITS_MI;
			}
		}
	
	float FinalM = Mrow[LB-1];
	float FinalD = Drow[LB];
	float FinalI = I1;
// FinalM = DPM[LA][LB]
// FinalD = DPD[LA][LB]
// FinalI = DPI[LA][LB]
	
	float Score = FinalM;
	byte State = 'M';
	if (FinalD > Score)
		{
		Score = FinalD;
		State = 'D';
		}
	if (FinalI > Score)
		{
		Score = FinalI;
		State = 'I';
		}

	TraceBackBitMem(Mem, LA, LB, State, Path);

	return Score;
	}
