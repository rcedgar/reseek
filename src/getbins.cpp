#include "myutils.h"
#include "dss.h"

void DSSParams::GetBins(FEATURE F, vector<float> &Bins)
	{
	Bins.clear();
	uint AlphaSize = GetAlphaSize(F);

#define BIN_T(Feat, Idx, t)	if (F == FEATURE_##Feat) Bins.push_back(float(t));

// Agree with DSS::ValueToInt_NormDens()
BIN_T(NormDens, 0, 0.241);
BIN_T(NormDens, 1, 0.3399);
BIN_T(NormDens, 2, 0.4115);
BIN_T(NormDens, 3, 0.4699);
BIN_T(NormDens, 4, 0.5204);
BIN_T(NormDens, 5, 0.5655);
BIN_T(NormDens, 6, 0.6065);
BIN_T(NormDens, 7, 0.6443);
BIN_T(NormDens, 8, 0.6803);
BIN_T(NormDens, 9, 0.715);
BIN_T(NormDens, 10, 0.7496);
BIN_T(NormDens, 11, 0.7854);
BIN_T(NormDens, 12, 0.8233);
BIN_T(NormDens, 13, 0.8655);
BIN_T(NormDens, 14, 0.917);

// Agree with DSS::ValueToInt_HelixDens()
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

// Agree with DSS::ValueToInt_StrandDens()
BIN_T(StrandDens, 0, 0.02212);
BIN_T(StrandDens, 1, 0.07567);
BIN_T(StrandDens, 2, 0.1134);
BIN_T(StrandDens, 3, 0.1394);
BIN_T(StrandDens, 4, 0.1605);
BIN_T(StrandDens, 5, 0.1796);
BIN_T(StrandDens, 6, 0.1982);
BIN_T(StrandDens, 7, 0.2172);
BIN_T(StrandDens, 8, 0.2378);
BIN_T(StrandDens, 9, 0.2615);
BIN_T(StrandDens, 10, 0.2893);
BIN_T(StrandDens, 11, 0.3227);
BIN_T(StrandDens, 12, 0.3627);
BIN_T(StrandDens, 13, 0.4111);
BIN_T(StrandDens, 14, 0.4778);

// Agree with DSS::ValueToInt_NENDist()
BIN_T(NENDist, 0, 4.417);
BIN_T(NENDist, 1, 4.647);
BIN_T(NENDist, 2, 4.841);
BIN_T(NENDist, 3, 5.052);
BIN_T(NENDist, 4, 5.286);
BIN_T(NENDist, 5, 5.589);
BIN_T(NENDist, 6, 6.055);
BIN_T(NENDist, 7, 6.536);
BIN_T(NENDist, 8, 7.007);
BIN_T(NENDist, 9, 7.485);
BIN_T(NENDist, 10, 7.999);
BIN_T(NENDist, 11, 8.559);
BIN_T(NENDist, 12, 9.166);
BIN_T(NENDist, 13, 9.873);
BIN_T(NENDist, 14, 11.18);

BIN_T(RENDist, 0, 6);
BIN_T(RENDist, 1, 7);
BIN_T(RENDist, 2, 8);
BIN_T(RENDist, 3, 9);
BIN_T(RENDist, 4, 10);
BIN_T(RENDist, 5, 11);
BIN_T(RENDist, 6, 12);
BIN_T(RENDist, 7, 13);
BIN_T(RENDist, 8, 14);
BIN_T(RENDist, 9, 15);
BIN_T(RENDist, 10, 16);
BIN_T(RENDist, 11, 17);
BIN_T(RENDist, 12, 18);
BIN_T(RENDist, 13, 19);
BIN_T(RENDist, 14, 20);

// Updated to agree with vartoint.cpp
BIN_T(DstNxtHlx, 0, 6);
BIN_T(DstNxtHlx, 1, 7);
BIN_T(DstNxtHlx, 2, 8);
BIN_T(DstNxtHlx, 3, 9);
BIN_T(DstNxtHlx, 4, 10);
BIN_T(DstNxtHlx, 5, 11);
BIN_T(DstNxtHlx, 6, 12);
BIN_T(DstNxtHlx, 7, 13);
BIN_T(DstNxtHlx, 8, 14);
BIN_T(DstNxtHlx, 9, 15);
BIN_T(DstNxtHlx, 10, 16);
BIN_T(DstNxtHlx, 11, 18);
BIN_T(DstNxtHlx, 12, 20);
BIN_T(DstNxtHlx, 13, 24);
BIN_T(DstNxtHlx, 14, 28);

// NOT IN vartoint.cpp
BIN_T(DstPrvHlx, 0, 6);
BIN_T(DstPrvHlx, 1, 7);
BIN_T(DstPrvHlx, 2, 8);
BIN_T(DstPrvHlx, 3, 9);
BIN_T(DstPrvHlx, 4, 10);
BIN_T(DstPrvHlx, 5, 11);
BIN_T(DstPrvHlx, 6, 12);
BIN_T(DstPrvHlx, 7, 13);
BIN_T(DstPrvHlx, 8, 14);
BIN_T(DstPrvHlx, 9, 15);
BIN_T(DstPrvHlx, 10, 16);
BIN_T(DstPrvHlx, 11, 18);
BIN_T(DstPrvHlx, 12, 20);
BIN_T(DstPrvHlx, 13, 24);
BIN_T(DstPrvHlx, 14, 28);

// Agrees with vartoint.cpp
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

// Agrees with vartoint.cpp
BIN_T(NX, 0, 20.65);
BIN_T(NX, 1, 23.54);
BIN_T(NX, 2, 25.62);
BIN_T(NX, 3, 27.43);
BIN_T(NX, 4, 29.14);
BIN_T(NX, 5, 30.76);
BIN_T(NX, 6, 32.3);
BIN_T(NX, 7, 33.78);
BIN_T(NX, 8, 35.22);
BIN_T(NX, 9, 36.61);
BIN_T(NX, 10, 37.96);
BIN_T(NX, 11, 39.34);
BIN_T(NX, 12, 40.77);
BIN_T(NX, 13, 42.39);
BIN_T(NX, 14, 44.47);

#undef BIN_T
	asserta(SIZE(Bins) + 1 == AlphaSize);
	}
