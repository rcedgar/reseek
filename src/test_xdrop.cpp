#include "myutils.h"
#include "xdpmem.h"

void SetBLOSUM62();
float GetBlosum62Score(char a, char b);
void MergeFwdBwd(uint LA, uint LB,
  uint FwdLoA, uint FwdLoB, const string &FwdPath,
  uint BwdHiA, uint BwdHiB, const string &BwdPath,
  uint &LoA, uint &LoB, uint &HiA, uint &HiB, string &Path);

const string *ptrA;
const string *ptrB;

static float SubFn(void *UserData, uint PosA, uint PosB)
	{
	asserta(PosA < SIZE(*ptrA));
	asserta(PosB < SIZE(*ptrB));
	char a = (*ptrA)[PosA];
	char b = (*ptrB)[PosB];
	float Score = GetBlosum62Score(a, b);
	return Score;
	}

static void LogAln(const string &A, const string &B, uint LoA, uint LoB,
  fn_SubstScore SubFn, float Open, float Ext, const string &Path)
	{
	if (Path.empty())
		return;
	uint PosA = LoA;
	uint PosB = LoB;
	uint ColCount = SIZE(Path);
	string RowA;
	string RowB;
	float Score = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			Score += SubFn(0, PosA, PosB);
			RowA += A[PosA++];
			RowB += B[PosB++];
			break;

		case 'D':
			if (Col == 0)
				Score += Open;
			else
				{
				if (Path[Col-1] == 'D')
					Score += Ext;
				else
					Score += Open;
				}
			RowA += A[PosA++];
			RowB += '-';
			break;

		case 'I':
			if (Col == 0)
				Score += Open;
			else
				{
				if (Path[Col-1] == 'I')
					Score += Ext;
				else
					Score += Open;
				}
			RowA += '-';
			RowB += B[PosB++];
			break;
			}
		}
	Log("\n");
	Log("%s\n", RowA.c_str());
	Log("%s\n", RowB.c_str());
	Log("Score %.3g\n", Score);
	}

static void Test(const string &A, const string &B)
	{
	XDPMem Mem;
	float Open = -3;
	float Ext = -1;
	float X = 8;

	ptrA = &A;
	ptrB = &B;
	const uint LA = SIZE(A);
	const uint LB = SIZE(B);

	Mx<float> SMx;
	SMx.Alloc("SMx", LA, LB);
	float **S = SMx.GetData();
	for (uint i = 0; i < LA; ++i)
		for (uint j = 0; j < LB; ++j)
			S[i][j] = SubFn(0, i, j);

	string SWPath;
	uint Loi, Loj, Leni, Lenj;
	Log("______________________________SWFast________________________\n");
	float SWScore =
	  SWFast(Mem, SMx, LA, LB, Open, Ext, Loi, Loj, Leni, Lenj, SWPath);
	uint LoA = Loi;
	uint LoB = Loj;
	ProgressLog("SW score = %.3g Path = %s\n", SWScore, SWPath.c_str());
	LogAln(A, B, LoA, LoB, SubFn, Open, Ext, SWPath);

	const uint ColCount = SIZE(SWPath);
	if (ColCount < 8)
		return;
	uint MidPosA = LoA;
	uint MidPosB = LoB;
	for (uint Col = 0; Col < ColCount/2; ++Col)
		{
		char c = SWPath[Col];
		if (c == 'M' || c == 'D') ++MidPosA;
		if (c == 'M' || c == 'I') ++MidPosB;
		}
	Log("Mid %u, %u\n", MidPosA, MidPosB);

	string FwdPath;
	Log("______________________________Fwd________________________\n");
	float FwdScore = XDropFwd(Mem, X, Open, Ext, SubFn, 0,
	  MidPosA+1, LA, MidPosB+1, LB, FwdPath);
	ProgressLog("FwdScore = %.3g Path = (%u,%u) %s\n",
	  FwdScore, MidPosA+1, MidPosB+1, FwdPath.c_str());
	LogAln(A, B, MidPosA, MidPosB, SubFn, Open, Ext, FwdPath);

	Log("______________________________Bwd________________________\n");
	string BwdPath;
	float BwdScore = XDropBwd(Mem, X, Open, Ext, SubFn, 0,
	  MidPosA, LA, MidPosB, LB, BwdPath);
	ProgressLog("BwdScore = %.3g (%u,%u) Path = %s\n",
	  BwdScore, MidPosA, MidPosB, BwdPath.c_str());
	uint LoLoA = MidPosA+1;
	uint LoLoB = MidPosB+1;
	uint ColCountB = SIZE(BwdPath);
	for (uint Col = 0; Col < ColCountB; ++Col)
		{
		char c = BwdPath[Col];
		if (c == 'M' || c == 'D')
			{
			asserta(LoLoA > 0);
			--LoLoA;
			}
		if (c == 'M' || c == 'I')
			{
			asserta(LoLoB > 0);
			--LoLoB;
			}
		}
	LogAln(A, B, LoLoA, LoLoB, SubFn, Open, Ext, BwdPath);
	char Mida = A[MidPosA];
	char Midb = B[MidPosB];
	float CombinedScore = FwdScore + BwdScore - GetBlosum62Score(Mida, Midb);
	string CombinedPath = BwdPath + FwdPath.substr(1);
	ProgressLog("FB score %.3g  %s\n", CombinedScore, CombinedPath.c_str());
	ProgressLog("SW score %.3g  %s\n", SWScore, SWPath.c_str());

	Log("______________________________Merged________________________\n");
	uint MergedLoA, MergedLoB;
	uint MergedHiA, MergedHiB;
	string MergedPath;
	MergeFwdBwd(LA, LB,
	  MidPosA+1, MidPosB+1, FwdPath,
	  MidPosA, MidPosB, BwdPath, MergedLoA, MergedLoB, MergedHiA, MergedHiB, MergedPath);
	Log("Merged A %u-%u, B %u-%u, Path %s\n",
	  MergedLoA, MergedLoB, MergedHiA, MergedHiB, MergedPath.c_str());
	LogAln(A, B, MergedLoA, MergedLoB, SubFn, Open, Ext, MergedPath);
	Log("====================================================================\n");
	}

void cmd_test_xdrop()
	{
	Test("DVLGYLRFLTKGERQANLNF",
		 "WVLGLRFLTKGERQANLNF");

	Test("DVLGYLRFLTERQANLNF",
		 "WVLGLRFLTKGERQANLNF");

	Test("DVLGYLRFLTKGERQANLNF",
		 "WVLGLINSRFLTKGERQANLNF");
	}
