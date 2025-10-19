#include "myutils.h"
#include "dss.h"

/***
Sequester Mu-related parameters and code here,
unrelated parameter tuning should not change Mu.

SEE ALSO
mumx_data.cpp
***/

char GetFeatureChar(byte Letter, uint AlphaSize);

static const uint s_MuFeatureCount = 3;
static const FEATURE s_MuFeatures[s_MuFeatureCount] = 
	{
	FEATURE_SS3,
	FEATURE_NENSS3,
	FEATURE_RENDist4
	};

static const uint s_MuAlphaSizes[3] = {3, 3, 4};;
static const uint s_MuAlphaSize = 36;
static float s_DefaultNENDist = 10.0;

static float RENDist_Ts[15] =
	{
#define BIN_T(Feat, Idx, t)	float(t),
BIN_T(RENDist, 0, 6)
BIN_T(RENDist, 1, 7)
BIN_T(RENDist, 2, 8)
BIN_T(RENDist, 3, 9)
BIN_T(RENDist, 4, 10)
BIN_T(RENDist, 5, 11)
BIN_T(RENDist, 6, 12)
BIN_T(RENDist, 7, 13)
BIN_T(RENDist, 8, 14)
BIN_T(RENDist, 9, 15)
BIN_T(RENDist, 10, 16)
BIN_T(RENDist, 11, 17)
BIN_T(RENDist, 12, 18)
BIN_T(RENDist, 13, 19)
BIN_T(RENDist, 14, 20)
#undef BIN_T
	};

static uint ValueToInt_RENDist_ForMu(float Value)
	{
	for (uint i = 0; i < 15; ++i)
		if (Value <= RENDist_Ts[i])
			return i;
	return 15;
	}

float DSS::GetFloat_RENDist_ForMu(uint Pos)
	{
	uint NEN = GetREN(Pos);
	if (NEN == UINT_MAX)
		{
		if (opt(force_undef))
			return FLT_MAX;
		return s_DefaultNENDist;
		}
	float d = m_Chain->GetDist(Pos, NEN);
	return d;
	}

uint DSS::Get_RENDist4(uint Pos)
	{
	float d = GetFloat_RENDist_ForMu(Pos);
	uint ND = ValueToInt_RENDist_ForMu(d);
	if (ND == 0)
		{
		if (opt(force_undef))
			return UINT_MAX;
		return UNDEFINED_ZERO_OVERLOAD;
		}
	asserta(ND < 16);
	return ND/4;
	}

uint DSS::Get_NENSS3(uint Pos)
	{
	SetSS();
	uint NEN = GetNEN(Pos);
	if (NEN == UINT_MAX)
		{
		if (opt(force_undef))
			return UINT_MAX;
		return 0;
		}
	char c = m_SS[NEN];
	switch (c)
		{
	case 'h': return 0;
	case 's': return 1;
	case 't': return 2;
	case '~': return 2;
		}
	asserta(false);
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
	asserta(false);
	return 0;
	}

void DSS::GetMuLetters(uint MuLetter, vector<uint> &Letters) const
	{
	Letters.clear();
	uint n = s_MuFeatureCount;
	if (MuLetter == UINT_MAX)
		{
		for (uint i = 0; i < n; ++i)
			Letters.push_back(UINT_MAX);
		return;
		}

	uint CL = MuLetter;
	uint m = 1;
	for (uint i = 0; i < n; ++i)
		{
		uint m = s_MuAlphaSizes[i];
		uint Letter = CL%m;
		Letters.push_back(Letter);
		CL /= m;
		}
	}

void DSS::GetMuLetters(vector<uint> &Letters)
	{
	Letters.clear();
	const uint L = GetSeqLength();
	Letters.reserve(L);
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		uint Letter = GetFeature(FEATURE_Mu, Pos);
		asserta(Letter < s_MuAlphaSize || Letter == UINT_MAX);
		Letters.push_back(Letter);
		}
	}

void DSS::GetMuKmers(const vector<byte> &Letters,
  vector<uint> &Kmers, const string &PatternStr)
	{
	Kmers.clear();
	const uint PatternLength = SIZE(PatternStr);
	const uint L = SIZE(Letters);
	Kmers.reserve(L);
	for (uint Pos = 0; Pos + PatternLength <= L; ++Pos)
		{
		uint Kmer = 0;
		for (uint j = 0; j < PatternLength; ++j)
			{
			if (PatternStr[j] == '1')
				{
				asserta(Pos + j < SIZE(Letters));
				uint Letter = Letters[Pos + j];
				assert(Letter < 36);
				Kmer = Kmer*36 + Letter;
				}
			}
		assert(Kmer != UINT_MAX);
		Kmers.push_back(Kmer);
		}
	}

uint DSS::Get_Mu(uint Pos)
	{
	uint MuLetter = 0;
	uint m = 1;
	for (uint i = 0; i < s_MuFeatureCount; ++i)
		{
		uint Letter = GetFeature(s_MuFeatures[i], Pos);
		if (Letter == UINT_MAX)
			return UINT_MAX;
		MuLetter = MuLetter + m*Letter;
		m *= s_MuAlphaSizes[i];
		}

	asserta(MuLetter < s_MuAlphaSize);
	return MuLetter;
	}

void cmd_mu_mapping()
	{
	DSS D;
	uint AS = s_MuAlphaSize;
	const uint N = s_MuFeatureCount;
	Log("Mu");
	for (uint i = 0; i < N; ++i)
		Log("\t%s", FeatureToStr(s_MuFeatures[i]));
	Log("\n");

	vector<uint> ASs;
	for (uint i = 0; i < N; ++i)
		ASs.push_back(DSS::GetAlphaSize(s_MuFeatures[i]));

	for (uint Letter = 0; Letter < AS; ++Letter)
		{
		vector<uint> Letters;
		D.GetMuLetters(Letter, Letters);
		asserta(SIZE(Letters) == N);
		char c = GetFeatureChar(Letter, AS);
		Log("%c", c);
		for (uint i = 0; i < N; ++i)
			{
			uint AS = ASs[i];
			uint Letter = Letters[i];
			asserta(Letter < AS);
			char c = GetFeatureChar(Letter, AS);
			Log("\t%c", c);
			}
		Log("\n");
		}
	}
