#include "myutils.h"
#include "mx.h"
#include "xdpmem.h"
#include "tracebit.h"
#include "pathinfo.h"
#include "swtrace.h"

struct RevData
	{
	uint LA;
	uint LB;
	void *ptrUserData = 0;
	ptr_fn_SubstScore ptrSubFn;
	};

float RevSubFn(void *UserData, uint RevPosA, uint RevPosB)
	{
	const RevData &RD = *(RevData *) UserData;
	assert(RevPosA < RD.LA);
	assert(RevPosB < RD.LB);
	uint PosA = RD.LA - RevPosA - 1;
	uint PosB = RD.LB - RevPosB - 1;
	void *ptrUserData = RD.ptrUserData;
	float Score = (*RD.ptrSubFn)(ptrUserData, PosA, PosB);
	return Score;
	}

float XDropBwd(XDPMem &Mem,
  float X, float Open, float Ext, 
  fn_SubstScore SubFn, void *UserData, 
  uint HiA, uint LA, uint HiB, uint LB,
  uint *ptrSegLoA, uint *ptrSegLoB,
  string &Path)
	{
	*ptrSegLoA = UINT_MAX;
	*ptrSegLoB = UINT_MAX;
	asserta(HiA < LA);
	asserta(HiB < LB);
	RevData RD;
	RD.LA = HiA+1;
	RD.LB = HiB+1;
	RD.ptrUserData = UserData;
	RD.ptrSubFn = SubFn;

	uint SegLoA, SegLoB;
	float Score = XDropFwd(Mem, X, Open, Ext, RevSubFn, &RD,
	  0, HiA+1, 0, HiB+1, &SegLoA, &SegLoB, Path);
	reverse(Path.begin(), Path.end());
	return Score;
	}
