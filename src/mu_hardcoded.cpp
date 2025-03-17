#include "myutils.h"
#include "dss.h"

static float DEFAULT_NEN_DIST = 10;

static uint ValueToInt_RENDist16(float Value)
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

uint DSS::Get_NENSS3(uint Pos)
	{
	SetSS();
	uint NEN = GetNEN(Pos);
	if (NEN == UINT_MAX)
		return 0;
	char c = m_SS[NEN];
	switch (c)
		{
	case 'h': return 0;
	case 's': return 1;
	case 't': return 2;
	case '~': return 2;
		}
	return 0;
	}

uint DSS::Get_SS3(uint Pos)
	{
	SetSS();
	char c = m_SS[Pos];
	switch (c)
		{
	case 'h': return 0;
	case 's': return 1;
	case 't': return 2;
	case '~': return 2;
		}
	return 0;
	}

//	SS3,
//	NENSS3,
//	RENDist4
//int m_NEN_W = 100;
//int m_NEN_w = 12;
uint DSS::Get_Mu(uint Pos)
	{
	//Log("Get_Mu(%u) ", Pos);//@@
	//uint NEN = GetNEN(Pos);
	//uint REN = GetREN(Pos);
	//if (NEN == UINT_MAX) Log(" NEN=*"); else Log(" NEN=%u", NEN);//@@
	//if (REN == UINT_MAX) Log(" REN=*"); else Log(" REN=%u", REN);//@@
	uint Letter_SS3 = Get_SS3(Pos);
	//Log(" SS3=%u", Letter_SS3);//@@
	assert(Letter_SS3 < 3);

	uint Letter_NENSS3 = Get_NENSS3(Pos);
	assert(Letter_NENSS3 < 3);
	//Log(" NENSS3=%u", Letter_NENSS3);//@@

	uint Letter_RENDist4 = 0;
	float RENDist = GetFloat_RENDist(Pos);
	if (RENDist == FLT_MAX)
		RENDist = DEFAULT_NEN_DIST;
	Letter_RENDist4 = ValueToInt_RENDist16(RENDist)/4;
	assert(Letter_RENDist4 < 4);
	//Log(" RENDist4=%u", Letter_RENDist4);//@@

	//uint MuLetter = Letter_SS3*3*4 + Letter_NENSS3*4 + Letter_RENDist4;
	uint MuLetter = Letter_SS3 + Letter_NENSS3*3 + Letter_RENDist4*3*3;
	assert(MuLetter < 36);
	//Log(" [%u]\n", MuLetter);//@@

	return MuLetter;
	}
