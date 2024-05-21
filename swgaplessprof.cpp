#include "myutils.h"
#include "mx.h"
#include "alpha.h"
#include "xdpmem.h"

extern float B62Mf[20][20];

// float GetBlosum62Score(char a, char b);
float SWGapless(Mx<float> &DPMx, const Mx<float> &SMx, uint LA, uint LB,
  uint &Loi, uint &Loj, uint &ColCount);

float SWFastGaplessProf(XDPMem &Mem, const float * const *ProfA, uint LA,
  const byte *B, uint LB, uint &Besti, uint &Bestj)
	{
	Mem.Alloc(LA+1, LB+1);

	Besti = UINT_MAX;
	Bestj = UINT_MAX;

	float *Mrow = Mem.GetDPRow1();

// Use Mrow[-1], so...
	Mrow[-1] = 0;
	for (uint j = 0; j <= LB; ++j)
		Mrow[j] = 0;
	
	float BestScore = 0.0f;

// Main loop
	float M0 = float (0);
	for (uint i = 0; i < LA; ++i)
		{
		const float *SMxRow = ProfA[i];
		for (uint j = 0; j < LB; ++j)
			{
			float SavedM0 = M0;

		// MATCH
			{
			float xM = M0;
			if (xM < 0.0f)
				xM = 0.0f;

			M0 = Mrow[j];
			byte bj = B[j];
			assert(bj < 36);
			xM += SMxRow[bj];
			if (xM > BestScore)
				{
				BestScore = xM;
				Besti = i;
				Bestj = j;
				}

			Mrow[j] = xM;
			}
			}
		
		M0 = 0;
		}
	
	return BestScore;
	}

static void MakeBlosumS(const string &A, const string &B,
  Mx<float> &MxS)
	{
	uint LA = SIZE(A);
	uint LB = SIZE(B);
	MxS.Alloc("BlosumS", LA, LB);
#if DEBUG
	MxS.Assign(FLT_MAX);
#endif
	float **S = MxS.GetData();
	for (uint i = 0; i < LA; ++i)
		{
		char a = A[i];
		byte ai = g_CharToLetterAmino[a];
		asserta(ai < 20);
		for (uint j = 0; j < LB; ++j)
			{
			char b = B[j];
			byte bi = g_CharToLetterAmino[b];
			asserta(bi < 20);
			S[i][j] = B62Mf[ai][bi];
			}
		}
	}

static void MakeProf(const string &A, vector<const float *> &Prof)
	{
	uint LA = SIZE(A);
	Prof.resize(LA);
	for (uint i = 0; i < LA; ++i)
		{
		char a = A[i];
		uint ai = g_CharToLetterAmino[a];
		asserta(ai < 20);
		Prof[i] = B62Mf[ai];
		}
	}

static void LogAln(const string &A, const string &B,
  uint Loi, uint Loj, uint ColCount)
	{
	const uint LA = SIZE(A);
	const uint LB = SIZE(B);
	asserta(Loi + ColCount <= LA);
	asserta(Loj + ColCount <= LB);
	Log("A ");
	for (uint Col = 0; Col < ColCount; ++Col)
		Log("%c", A[Loi+Col]);
	Log(" %u .. %u\n", Loi + 1, Loi + ColCount);

	Log("B ");
	for (uint Col = 0; Col < ColCount; ++Col)
		Log("%c", B[Loj+Col]);
	Log(" %u .. %u\n", Loj + 1, Loj + ColCount);
	}

static void Test2(const string &A, const string &B)
	{
	uint LA = SIZE(A);
	uint LB = SIZE(B);

	Mx<float> SMx;
	Mx<float> DPMx;
	MakeBlosumS(A, B, SMx);

	vector<const float *> ProfA;
	MakeProf(A, ProfA);

	vector<byte> Bi;
	for (uint i = 0; i < LB; ++i)
		{
		char c = B[i];
		byte ci = g_CharToLetterAmino[c];
		asserta(ci < 20);
		Bi.push_back(ci);
		}

	uint Loi, Loj, ColCount;
	float Score = SWGapless(DPMx, SMx, LA, LB, Loi, Loj, ColCount);
	Log("\n");
	Log("A=%s\n", A.c_str());
	Log("B=%s\n", B.c_str());
	Log("Score %.1f, %u, %u, %u\n", Score, Loi, Loj, ColCount);
	if (Score > 0)
		LogAln(A, B, Loi, Loj, ColCount);

	XDPMem Mem;
	uint Besti, Bestj;
	float Score2 = SWFastGaplessProf(Mem, ProfA.data(), LA, Bi.data(), LB, Besti, Bestj);
	Log("Score2 %.1f, %u, %u\n", Score2, Besti, Bestj);
	asserta(feq(Score, Score2));
	}

void cmd_testsw()
	{
	Test2("LQNGSEQVENCE", "LQNGSEQVENCE");
	Test2("QNGSEQVENCE", "LQNGSEQVENCE");
	Test2("LQNGSEQVENCE", "QNGSEQVENCE");
	Test2("LQNGSEQVENC", "QNGSEQVENCE");
	Test2("SEQVENCE", "QVE");
	Test2("QVE", "SEQVENCE");
	}
