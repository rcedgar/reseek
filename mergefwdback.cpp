#include "myutils.h"
#include "pdbchain.h"

void GetPathCounts(const string &Path, uint &M, uint &D, uint &I)
	{
	M = D = I = 0;
	const uint ColCount = SIZE(Path);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		switch (Path[Col])
			{
		case 'M': ++M; continue;
		case 'D': ++D; continue;
		case 'I': ++I; continue;
		default: asserta(false);
			}
		}
	}

void MergeFwdBwd(uint LA, uint LB,
  uint FwdLoA, uint FwdLoB, const string &FwdPath,
  uint BwdHiA, uint BwdHiB, const string &BwdPath,
  uint &LoA, uint &LoB, uint &HiA, uint &HiB, string &Path)
	{
	asserta(FwdPath != "" || BwdPath != "");
	asserta(FwdLoA == BwdHiA + 1);
	asserta(FwdLoB == BwdHiB + 1);

	if (FwdPath == "")
		{
		HiA = BwdHiA;
		HiB = BwdHiB;
		}
	else
		{
		uint FwdM, FwdD, FwdI;
		GetPathCounts(FwdPath, FwdM, FwdD, FwdI);
		uint FwdColsA = FwdM + FwdD;
		uint FwdColsB = FwdM + FwdI;
		HiA = FwdLoA + FwdColsA - 1;
		HiB = FwdLoB + FwdColsB - 1;
		asserta(HiA < LA);
		asserta(HiB < LB);
		}

	if (BwdPath == "")
		{
		LoA = FwdLoA;
		LoB = FwdLoB;
		}
	else
		{
		uint BwdM, BwdD, BwdI;
		GetPathCounts(BwdPath, BwdM, BwdD, BwdI);

		uint BwdColsA = BwdM + BwdD;
		uint BwdColsB = BwdM + BwdI;
		asserta(BwdHiA + 1 >= BwdColsA);
		asserta(BwdHiB + 1 >= BwdColsB);
		LoA = BwdHiA + 1 - BwdColsA;
		LoB = BwdHiB + 1 - BwdColsB;
		}
	Path = BwdPath + FwdPath;
	}
