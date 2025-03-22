#include "myutils.h"
#include "dss.h"

#if 0

#define BIN_T(Feat, Idx, t)	if (Value < t) return Idx;

uint DSS::ValueToInt_NENDist(float Value) const
	{
	if (Value < 4.417) return 0;
	if (Value < 4.647) return 1;
	if (Value < 4.841) return 2;
	if (Value < 5.052) return 3;
	if (Value < 5.286) return 4;
	if (Value < 5.589) return 5;
	if (Value < 6.055) return 6;
	if (Value < 6.536) return 7;
	if (Value < 7.007) return 8;
	if (Value < 7.485) return 9;
	if (Value < 7.999) return 10;
	if (Value < 8.559) return 11;
	if (Value < 9.166) return 12;
	if (Value < 9.873) return 13;
	if (Value < 11.18) return 14;
	return 15;
	}

uint DSS::ValueToInt_RENDist(float Value) const
	{
	if (Value < 6) return 0;
	if (Value < 7) return 1;
	if (Value < 8) return 2;
	if (Value < 9) return 3;
	if (Value < 10) return 4;
	if (Value < 11) return 5;
	if (Value < 12) return 6;
	if (Value < 13) return 7;
	if (Value < 14) return 8;
	if (Value < 15) return 9;
	if (Value < 16) return 10;
	if (Value < 17) return 11;
	if (Value < 18) return 12;
	if (Value < 19) return 13;
	if (Value < 20) return 14;
	return 15;
	}

uint DSS::ValueToInt_DstNxtHlx(float Value) const
	{
	if (Value < 6) return 0;
	if (Value < 7) return 1;
	if (Value < 8) return 2;
	if (Value < 9) return 3;
	if (Value < 10) return 4;
	if (Value < 11) return 5;
	if (Value < 12) return 6;
	if (Value < 13) return 7;
	if (Value < 14) return 8;
	if (Value < 15) return 9;
	if (Value < 16) return 10;
	if (Value < 18) return 11;
	if (Value < 20) return 12;
	if (Value < 24) return 13;
	if (Value < 28) return 14;
	return 15;
	}

uint DSS::ValueToInt_StrandDens(float Value) const
	{
	if (Value < 0.02212) return 0;
	if (Value < 0.07567) return 1;
	if (Value < 0.1134) return 2;
	if (Value < 0.1394) return 3;
	if (Value < 0.1605) return 4;
	if (Value < 0.1796) return 5;
	if (Value < 0.1982) return 6;
	if (Value < 0.2172) return 7;
	if (Value < 0.2378) return 8;
	if (Value < 0.2615) return 9;
	if (Value < 0.2893) return 10;
	if (Value < 0.3227) return 11;
	if (Value < 0.3627) return 12;
	if (Value < 0.4111) return 13;
	if (Value < 0.4778) return 14;
	return 15;
	}

uint DSS::ValueToInt_NormDens(float Value) const
	{
	if (Value < 0.241) return 0;
	if (Value < 0.3399) return 1;
	if (Value < 0.4115) return 2;
	if (Value < 0.4699) return 3;
	if (Value < 0.5204) return 4;
	if (Value < 0.5655) return 5;
	if (Value < 0.6065) return 6;
	if (Value < 0.6443) return 7;
	if (Value < 0.6803) return 8;
	if (Value < 0.715) return 9;
	if (Value < 0.7496) return 10;
	if (Value < 0.7854) return 11;
	if (Value < 0.8233) return 12;
	if (Value < 0.8655) return 13;
	if (Value < 0.917) return 14;
	return 15;
	}

uint DSS::ValueToInt_HelixDens(float Value) const
	{
BIN_T(HelixDens, 0, 0.03015);
BIN_T(HelixDens, 1, 0.06112);
BIN_T(HelixDens, 2, 0.1127);
BIN_T(HelixDens, 3, 0.1683);
BIN_T(HelixDens, 4, 0.2115);
BIN_T(HelixDens, 5, 0.2455);
BIN_T(HelixDens, 6, 0.275);
BIN_T(HelixDens, 7, 0.3033);
BIN_T(HelixDens, 8, 0.3309);
BIN_T(HelixDens, 9, 0.3589);
BIN_T(HelixDens, 10, 0.3885);
BIN_T(HelixDens, 11, 0.4227);
BIN_T(HelixDens, 12, 0.4647);
BIN_T(HelixDens, 13, 0.5258);
BIN_T(HelixDens, 14, 0.6343);
	return 15;
	}

uint DSS::ValueToInt_PMDist(float Value) const
	{
BIN_T(PMDist, 0, 9.994);
BIN_T(PMDist, 1, 12.06);
BIN_T(PMDist, 2, 13.65);
BIN_T(PMDist, 3, 14.98);
BIN_T(PMDist, 4, 16.3);
BIN_T(PMDist, 5, 17.57);
BIN_T(PMDist, 6, 18.82);
BIN_T(PMDist, 7, 20.06);
BIN_T(PMDist, 8, 21.33);
BIN_T(PMDist, 9, 22.64);
BIN_T(PMDist, 10, 23.93);
BIN_T(PMDist, 11, 24.86);
BIN_T(PMDist, 12, 26.38);
BIN_T(PMDist, 13, 28.84);
BIN_T(PMDist, 14, 32.77);
	return 15;
	}

uint DSS::ValueToInt_DstPrvHlx(float Value) const
	{
BIN_T(DstNxtHlx, 0, 0);
BIN_T(DstNxtHlx, 1, 6);
BIN_T(DstNxtHlx, 2, 7);
BIN_T(DstNxtHlx, 3, 8);
BIN_T(DstNxtHlx, 4, 9);
BIN_T(DstNxtHlx, 5, 10.81);
BIN_T(DstNxtHlx, 6, 12.59);
BIN_T(DstNxtHlx, 7, 14.01);
BIN_T(DstNxtHlx, 8, 15.25);
BIN_T(DstNxtHlx, 9, 16.62);
BIN_T(DstNxtHlx, 10, 18.21);
BIN_T(DstNxtHlx, 11, 19.98);
BIN_T(DstNxtHlx, 12, 22);
BIN_T(DstNxtHlx, 13, 24.6);
BIN_T(DstNxtHlx, 14, 28.82);
	return 15;
	}

uint DSS::ValueToInt_PlusNENDist(float Value) const
	{
BIN_T(PlusNENDist, 0, 4.633);
BIN_T(PlusNENDist, 1, 4.937);
BIN_T(PlusNENDist, 2, 5.292);
BIN_T(PlusNENDist, 3, 5.924);
BIN_T(PlusNENDist, 4, 6.701);
BIN_T(PlusNENDist, 5, 7.436);
BIN_T(PlusNENDist, 6, 8.17);
BIN_T(PlusNENDist, 7, 8.861);
BIN_T(PlusNENDist, 8, 9.48);
BIN_T(PlusNENDist, 9, 10);
BIN_T(PlusNENDist, 10, 10.11);
BIN_T(PlusNENDist, 11, 11.05);
BIN_T(PlusNENDist, 12, 12.47);
BIN_T(PlusNENDist, 13, 14.69);
BIN_T(PlusNENDist, 14, 18.8);
	return 15;
	}

uint DSS::ValueToInt_MinusNENDist(float Value) const
	{
BIN_T(MinusNENDist, 0, 4.634);
BIN_T(MinusNENDist, 1, 4.92);
BIN_T(MinusNENDist, 2, 5.244);
BIN_T(MinusNENDist, 3, 5.788);
BIN_T(MinusNENDist, 4, 6.613);
BIN_T(MinusNENDist, 5, 7.348);
BIN_T(MinusNENDist, 6, 8.027);
BIN_T(MinusNENDist, 7, 8.698);
BIN_T(MinusNENDist, 8, 9.33);
BIN_T(MinusNENDist, 9, 9.964);
BIN_T(MinusNENDist, 10, 10);
BIN_T(MinusNENDist, 11, 10.76);
BIN_T(MinusNENDist, 12, 12.02);
BIN_T(MinusNENDist, 13, 14.14);
BIN_T(MinusNENDist, 14, 18.33);
	return 15;
	}

uint DSS::ValueToInt_DiffNENDist(float Value) const
	{
BIN_T(DiffNENDist, 0, -10);
BIN_T(DiffNENDist, 1, -6.277);
BIN_T(DiffNENDist, 2, -4.342);
BIN_T(DiffNENDist, 3, -2.889);
BIN_T(DiffNENDist, 4, -1.519);
BIN_T(DiffNENDist, 5, -0.4359);
BIN_T(DiffNENDist, 6, 0);
BIN_T(DiffNENDist, 7, 0.1911);
BIN_T(DiffNENDist, 8, 1.121);
BIN_T(DiffNENDist, 9, 2.368);
BIN_T(DiffNENDist, 10, 3.693);
BIN_T(DiffNENDist, 11, 4.978);
BIN_T(DiffNENDist, 12, 7.068);
BIN_T(DiffNENDist, 13, 10);
BIN_T(DiffNENDist, 14, 10.28);
	return 15;
	}

#endif // 0