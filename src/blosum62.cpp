#include "myutils.h"
#include "mx.h"
#include "alpha.h"
#include <mutex>

Mx<float> g_SubstMxf;
float **g_SubstMx;

/***
This alphabet is the one used by BLAST.
It includes wildcards (BZX) and stop codon (*).
This substitution matrix is used when alignments are constructed,
see uses of g_SubstMx.
The alphabet for word-counting is different, it has only the 20
standard amino acids. 
Keeping the word-counting alphabet smaller saves memory for
word index tables.
***/
static char Alphabet[] = "*ACBEDGFIHKMLNQPSRTWVYXZ";

// Ye olde BLOSUM62 as used by NCBI BLAST (1/2-bit units)
static float BLOSUM62[24][24] =
	{
//	    *    A    C    B    E    D    G    F    I    H    K    M    L    N    Q    P    S    R    T    W    V    Y    X    Z   
	{   1,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  },  // *
	{  -4,   4,   0,  -2,  -1,  -2,   0,  -2,  -1,  -2,  -1,  -1,  -1,  -2,  -1,  -1,   1,  -1,   0,  -3,   0,  -2,   0,  -1,  },  // A
	{  -4,   0,   9,  -3,  -4,  -3,  -3,  -2,  -1,  -3,  -3,  -1,  -1,  -3,  -3,  -3,  -1,  -3,  -1,  -2,  -1,  -2,  -2,  -3,  },  // C
	{  -4,  -2,  -3,   4,   1,   4,  -1,  -3,  -3,   0,   0,  -3,  -4,   3,   0,  -2,   0,  -1,  -1,  -4,  -3,  -3,  -1,   1,  },  // B
	{  -4,  -1,  -4,   1,   5,   2,  -2,  -3,  -3,   0,   1,  -2,  -3,   0,   2,  -1,   0,   0,  -1,  -3,  -2,  -2,  -1,   4,  },  // E
	{  -4,  -2,  -3,   4,   2,   6,  -1,  -3,  -3,  -1,  -1,  -3,  -4,   1,   0,  -1,   0,  -2,  -1,  -4,  -3,  -3,  -1,   1,  },  // D
	{  -4,   0,  -3,  -1,  -2,  -1,   6,  -3,  -4,  -2,  -2,  -3,  -4,   0,  -2,  -2,   0,  -2,  -2,  -2,  -3,  -3,  -1,  -2,  },  // G
	{  -4,  -2,  -2,  -3,  -3,  -3,  -3,   6,   0,  -1,  -3,   0,   0,  -3,  -3,  -4,  -2,  -3,  -2,   1,  -1,   3,  -1,  -3,  },  // F
	{  -4,  -1,  -1,  -3,  -3,  -3,  -4,   0,   4,  -3,  -3,   1,   2,  -3,  -3,  -3,  -2,  -3,  -1,  -3,   3,  -1,  -1,  -3,  },  // I
	{  -4,  -2,  -3,   0,   0,  -1,  -2,  -1,  -3,   8,  -1,  -2,  -3,   1,   0,  -2,  -1,   0,  -2,  -2,  -3,   2,  -1,   0,  },  // H
	{  -4,  -1,  -3,   0,   1,  -1,  -2,  -3,  -3,  -1,   5,  -1,  -2,   0,   1,  -1,   0,   2,  -1,  -3,  -2,  -2,  -1,   1,  },  // K
	{  -4,  -1,  -1,  -3,  -2,  -3,  -3,   0,   1,  -2,  -1,   5,   2,  -2,   0,  -2,  -1,  -1,  -1,  -1,   1,  -1,  -1,  -1,  },  // M
	{  -4,  -1,  -1,  -4,  -3,  -4,  -4,   0,   2,  -3,  -2,   2,   4,  -3,  -2,  -3,  -2,  -2,  -1,  -2,   1,  -1,  -1,  -3,  },  // L
	{  -4,  -2,  -3,   3,   0,   1,   0,  -3,  -3,   1,   0,  -2,  -3,   6,   0,  -2,   1,   0,   0,  -4,  -3,  -2,  -1,   0,  },  // N
	{  -4,  -1,  -3,   0,   2,   0,  -2,  -3,  -3,   0,   1,   0,  -2,   0,   5,  -1,   0,   1,  -1,  -2,  -2,  -1,  -1,   3,  },  // Q
	{  -4,  -1,  -3,  -2,  -1,  -1,  -2,  -4,  -3,  -2,  -1,  -2,  -3,  -2,  -1,   7,  -1,  -2,  -1,  -4,  -2,  -3,  -2,  -1,  },  // P
	{  -4,   1,  -1,   0,   0,   0,   0,  -2,  -2,  -1,   0,  -1,  -2,   1,   0,  -1,   4,  -1,   1,  -3,  -2,  -2,   0,   0,  },  // S
	{  -4,  -1,  -3,  -1,   0,  -2,  -2,  -3,  -3,   0,   2,  -1,  -2,   0,   1,  -2,  -1,   5,  -1,  -3,  -3,  -2,  -1,   0,  },  // R
	{  -4,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -2,  -1,  -1,  -1,   0,  -1,  -1,   1,  -1,   5,  -2,   0,  -2,   0,  -1,  },  // T
	{  -4,  -3,  -2,  -4,  -3,  -4,  -2,   1,  -3,  -2,  -3,  -1,  -2,  -4,  -2,  -4,  -3,  -3,  -2,  11,  -3,   2,  -2,  -3,  },  // W
	{  -4,   0,  -1,  -3,  -2,  -3,  -3,  -1,   3,  -3,  -2,   1,   1,  -3,  -2,  -2,  -2,  -3,   0,  -3,   4,  -1,  -1,  -2,  },  // V
	{  -4,  -2,  -2,  -3,  -2,  -3,  -3,   3,  -1,   2,  -2,  -1,  -1,  -2,  -1,  -3,  -2,  -2,  -2,   2,  -1,   7,  -1,  -2,  },  // Y
	{  -4,   0,  -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -2,   0,  -1,   0,  -2,  -1,  -1,  -1,  -1,  },  // X
	{  -4,  -1,  -3,   1,   4,   1,  -2,  -3,  -3,   0,   1,  -1,  -3,   0,   3,  -1,   0,   0,  -1,  -3,  -2,  -2,  -1,   4,  },  // Z
	};

void SetBLOSUM62Mx(Mx<float> &Sf)
	{
	unsigned N = unsigned(strlen(Alphabet));

	Sf.Alloc(256, 256, __FILE__, __LINE__);
	Sf.Init(0);
	float **Data = Sf.GetData();

	for (unsigned i = 0; i < N; ++i)
		{
		for (unsigned j = 0; j < N; ++j)
			{
			float v = BLOSUM62[i][j];
			
			byte ui = (byte) toupper(Alphabet[i]);
			byte uj = (byte) toupper(Alphabet[j]);
			byte li = (byte) tolower(ui);
			byte lj = (byte) tolower(uj);
			ui = (byte) toupper(ui);
			uj = (byte) toupper(uj);

			Data[ui][uj] = v;
			Data[uj][ui] = v;

			Data[ui][lj] = v;
			Data[uj][li] = v;

			Data[li][uj] = v;
			Data[lj][ui] = v;

			Data[li][lj] = v;
			Data[lj][li] = v;
			}
		}
	}

void SetBLOSUM62()
	{
	static mutex Lock;
	Lock.lock();
	if (g_SubstMx == 0)
		{
		SetBLOSUM62Mx(g_SubstMxf);
		g_SubstMx = g_SubstMxf.GetData();
		}
	Lock.unlock();
	}

float GetBlosum62Score(char a, char b)
	{
	static bool InitDone = false;
	if (!InitDone)
		{
		SetBLOSUM62();
		InitDone = true;
		}

	float Score = g_SubstMx[(byte) a][(byte) b];
	return Score;
	}

float GetBlosum62PathScore(
  const string &A, uint LoA,
  const string &B, uint LoB,
  const string &Path)
	{
	const float Open = -10;
	const float Ext = -1;
	float Sum = 0;
	uint LA = SIZE(A);
	uint LB = SIZE(B);
	uint PosA = LoA;
	uint PosB = LoB;
	const uint ColCount = SIZE(Path);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			{
			char a = A[PosA];
			char b = B[PosB];
			Sum += GetBlosum62Score(a, b);
			++PosA;
			++PosB;
			break;
			}

		case 'D':
			if (Col != 0 && Path[Col-1] == 'D')
				Sum += Ext;
			else
				Sum += Open;
			++PosA;
			break;

		case 'I':
			if (Col != 0 && Path[Col-1] == 'I')
				Sum += Ext;
			else
				Sum += Open;
			++PosB;
			break;

		default:
			asserta(false);
			}
		}
	asserta(PosA <= SIZE(A));
	asserta(PosB <= SIZE(B));
	return Sum;
	}

float GetBlastpEvalue(float Score, uint QL, float DBSize)
	{
	const float Lambda = 0.267f;
	const float K = 0.0410f;
	static const float Log2 = logf(2.0f);
	float BitScore = (Score*Lambda - K)/Log2;
	float NM = float(QL)*float(DBSize);
	float Evalue = NM/powf(2, BitScore);
	return Evalue;
	}
