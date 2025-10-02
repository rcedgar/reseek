#include "myutils.h"
#include "dss.h"

static uint ValueToInt_RENDist4(float Value)
	{
	if (Value < 9.194) return 0;
	if (Value < 11.74) return 1;
	if (Value < 16.66) return 2;
	return 3;
	}

uint DSS::Get_RENDist4(uint Pos)
	{
	float RENDist = GetFloat_RENDist(Pos);
	uint Letter = ValueToInt_RENDist4(RENDist);
	return Letter;
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

// -train_feature ../traindata/train.fa2 -db ../data/scop40.bca -feature NENSS3 -alpha_size 3 -wildcard no -log train_feature.log 	uint Letter_NENSS3 = Get_NENSS3(Pos);
//  0 = -0.143 << best default
//  1 = -0.344
//  2 = -0.164
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

// -train_feature ../traindata/train.fa2 -db ../data/scop40.bca -feature SS3 -alpha_size 3 -wildcard no -log train_feature.log 
// SS3(3) expected scores
//  0 = -1.29
//  1 = -0.703
//  2 = -0.628 << best default
	return 2;
	}

//	SS3,
//	NENSS3,
//	RENDist4
//int m_NEN_W = 100;
//int m_NEN_w = 12;
uint DSS::Get_Mu(uint Pos)
	{
	uint Letter_SS3 = Get_SS3(Pos);
	assert(Letter_SS3 < 3);
 
	uint Letter_NENSS3 = Get_NENSS3(Pos);
	assert(Letter_NENSS3 < 3);

// -train_feature ../traindata/train.fa2 -db ../data/scop40.bca -feature RENDist -alpha_size 4 -wildcard no -log train_feature.log 
// RENDist(4) expected scores
//  0 = -0.413
//  1 = -0.995
//  2 = -0.418
//  3 = -0.301 << best default
	uint Letter_RENDist4 = 3;
	float RENDist = GetFloat_RENDist(Pos);
	if (RENDist != FLT_MAX)
		Letter_RENDist4 = ValueToInt_RENDist4(RENDist);
	assert(Letter_RENDist4 < 4);

	uint MuLetter = Letter_SS3 + Letter_NENSS3*3 + Letter_RENDist4*3*3;
	assert(MuLetter < 36);
	return MuLetter;
	}
