#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "alpha.h"
#include "dss.h"

uint GetPatternOnes(const string &Str)
	{
	uint n = 0;
	for (uint i = 0; i < SIZE(Str); ++i)
		if (Str[i] == '1')
			++n;
	return n;
	}

void DSS::SetSS()
	{
	if (m_SS.empty())
		m_Chain->GetSS(m_SS);
	}

uint DSS::Get_SS(uint Pos)
	{
	SetSS();
	char c = m_SS[Pos];
	uint Letter = SSCharToInt(c);
	return Letter;
	}

uint DSS::Get_RENSS3(uint Pos)
	{
	SetSS();
	uint NEN = GetREN(Pos);
	if (NEN == UINT_MAX)
		return UNDEFINED_ZERO_OVERLOAD;
	char c = m_SS[NEN];
	switch (c)
		{
	case 'h': return 0;
	case 's': return 1;
	case 't': return 2;
	case '~': return 2;
		}
	return UNDEFINED_ZERO_OVERLOAD;
	}

void DSS::GetSSEs(uint MinLength, vector<uint> &Los,
  vector<uint> &Lengths, vector<char> &cs) const
	{
	Los.clear();
	Lengths.clear();
	cs.clear();
	const uint L = GetSeqLength();
	asserta(SIZE(m_SS) == L);
	char currc = m_SS[0];
	uint StartPos = 0;
	uint RunLength = 1;
	for (uint Pos = 1; Pos <= L; ++Pos)
		{
		char ss = (Pos == L ? 0 : m_SS[Pos]);
		if (m_SS[Pos] == currc)
			++RunLength;
		else
			{
			if (RunLength >= MinLength)
				{
				if (currc == 'h' || currc == 's')
					{
					Los.push_back(StartPos);
					Lengths.push_back(RunLength);
					cs.push_back(currc);
					}
				}
			currc = m_SS[Pos];
			StartPos = Pos;
			RunLength = 1;
			}
		}
	}

uint DSS::SSCharToInt(char c)
	{
	switch (c)
		{
	case 'h': return 0;
	case 's': return 1;
	case 't': return 2;
	case '~': return 3;
		}
	asserta(false);
	return UINT_MAX;
	}

uint DSS::SSCharToInt3(char c)
	{
	switch (c)
		{
	case 'h': return 0;
	case 's': return 1;
	case 't': return 2;
	case '~': return 2;
		}
	asserta(false);
	return UINT_MAX;
	}

void DSS::SetSSEs()
	{
	if (!m_SSE_Mids.empty())
		return;
	asserta(m_SSE_cs.empty());
	SetSS();
	vector<uint> Los;
	vector<uint> Lengths;
	GetSSEs(DSSParams::m_SSE_MinLength, Los, Lengths, m_SSE_cs);
	const uint N = SIZE(Los);
	asserta(SIZE(Lengths) == N);
	asserta(SIZE(m_SSE_cs) == N);
	for (uint i = 0; i < N; ++i)
		{
		uint Mid = Los[i] + Lengths[i]/2;
		m_SSE_Mids.push_back(Mid);
		}
	}

float DSS::GetFloat_NormDens(uint Pos)
	{
	SetDensity_ScaledValues();
	asserta(Pos < SIZE(m_Density_ScaledValues));
	float Valued = m_Density_ScaledValues[Pos];
	return Valued;
	}

float DSS::GetFloat_HelixDens(uint Pos)
	{
	return (float) GetSSDensity(Pos, 'h');
	}

float DSS::GetFloat_StrandDens(uint Pos)
	{
	return (float) GetSSDensity(Pos, 's');
	}

void DSS::SetDensity_ScaledValues()
	{
	if (!m_Density_ScaledValues.empty())
		return;
	m_Density_ScaledValues.clear();
	const uint L = GetSeqLength();
	vector<float> Values;
	vector<float> Valuesf;
	m_Density_ScaledValues.reserve(L);
	Values.reserve(L);
	float MinValue = 999;
	float MaxValue = 0;
	float MinValuef = 999;
	float MaxValuef = 0;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		float D = GetDensity(Pos);
		Values.push_back(D);
		if (D == FLT_MAX)
			Valuesf.push_back(FLT_MAX);
		else
			Valuesf.push_back(float(D));
		if (D != FLT_MAX)
			{
			MinValue = min(MinValue, D);
			MaxValue = max(MaxValue, D);
			MinValuef = min(MinValuef, float(D));
			}
		}

	float Range = (MaxValue - MinValue);
	float Rangef = (MaxValuef - MinValuef);
	if (Range < 1)
		Range = 1;
	if (Rangef < 1)
		Rangef = 1;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		float Value = Values[Pos];
		float Valuef = Valuesf[Pos];
		if (Value == FLT_MAX)
			{
			asserta(Valuef == FLT_MAX);
			m_Density_ScaledValues.push_back(FLT_MAX);
			continue;
			}
		float ScaledValue = (Value - MinValue)/Range;
		asserta(ScaledValue >= 0 && ScaledValue <= 1);
		m_Density_ScaledValues.push_back(ScaledValue);
		}
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		float Value = Values[Pos];
		float Valuef = Valuesf[Pos];
		if (Value == FLT_MAX)
			asserta(Valuef == FLT_MAX);
		else
			asserta(feq(Value, Valuef));
		}
	}

float DSS::GetDensity(uint Pos) const
	{
	const PDBChain &Chain = *m_Chain;
	const uint L = SIZE(Chain.m_Seq);
	if (Pos == 0 || Pos+1 >= L)
		return FLT_MAX;

	//vector<float> PtCA;
	//Chain.GetPt(Pos, PtCA);

	int iLo = int(Pos) - DSSParams::m_Density_W;
	if (iLo < 0)
		iLo = 0;
	int iHi = int(Pos) + DSSParams::m_Density_W;
	if (iHi >= int(L))
		iHi = int(L)-1;
	float D = 0;
	//vector<float> Pt2;
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + DSSParams::m_Density_w >= Pos && Pos2 <= Pos + DSSParams::m_Density_w)
			continue;
		float Dist = Chain.GetDist(Pos, Pos2);
		float DistFactor = exp(-Dist/DSSParams::m_Density_Radius);
		D += DistFactor;
		}
	return D;
	}

void DSS::Get_NU_ND(uint Pos, float &NU, float &ND) const
	{
	NU = 0;
	ND = 0;
	const PDBChain &Chain = *m_Chain;
	const uint L = SIZE(Chain.m_Seq);
	if (Pos == 0 || Pos+1 >= L)
		{
		NU = FLT_MAX;
		ND = FLT_MAX;
		return;
		}

	vector<float> PtPrevCA;
	vector<float> PtNextCA;
	vector<float> PtCA;
	vector<float> PtCB;
	Chain.GetPt(Pos-1, PtPrevCA);
	Chain.GetPt(Pos, PtCA);
	Chain.GetPt(Pos+1, PtNextCA);

	vector<float> d1;
	vector<float> d2;
	Sub_Vecs(PtCA, PtPrevCA, d1);
	Sub_Vecs(PtCA, PtNextCA, d2);

	vector<float> VecPAB;
	Add_Vecs(d1, d2, VecPAB);
	NormalizeVec(VecPAB);

	vector<float> Pt2;
	vector<float> Vec12;
	//const int W = 50;
	//const float RADIUS = 20.0;
	int iLo = int(Pos) - DSSParams::m_NUDX_W;
	if (iLo < 0)
		iLo = 0;
	int iHi = int(Pos) + DSSParams::m_NUDX_W;
	if (iHi >= int(L))
		iHi = int(L)-1;
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + 3 >= Pos && Pos2 <= Pos + 3)
			continue;
		float Dist = Chain.GetDist(Pos, Pos2);
		float DistFactor = exp(-Dist/DSSParams::m_NU_ND_Radius);
		Chain.GetPt(Pos2, Pt2);
		Sub_Vecs(Pt2, PtCA, Vec12);
		float Theta = GetTheta_Vecs(VecPAB, Vec12);
		float Deg = degrees(Theta);
		if (Deg < 90)
			NU += DistFactor;
		else
			ND += DistFactor;
		}
	}

void DSS::Set_NU_ND_Vecs()
	{
	if (!m_NUs.empty())
		return;
	const uint L = GetSeqLength();
	m_NUs.reserve(L);
	m_NDs.reserve(L);
	m_NXs.reserve(L);
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		float NU, ND;
		Get_NU_ND(Pos, NU, ND);
		m_NUs.push_back(NU);
		m_NDs.push_back(ND);
		if (NU == FLT_MAX || ND == FLT_MAX)
			m_NXs.push_back(FLT_MAX);
		else	
			m_NXs.push_back(NU+ND);
		}
	}

float DSS::GetFloat_NX(uint Pos)
	{
	Set_NU_ND_Vecs();
	return (float) m_NXs[Pos];
	}

//float DSS::GetFloat_NU(uint Pos)
//	{
//	Set_NU_ND_Vecs();
//	return m_NUs[Pos];
//	}
//
//float DSS::GetFloat_ND(uint Pos)
//	{
//	Set_NU_ND_Vecs();
//	return m_NDs[Pos];
//	}

float DSS::GetSSDensity(uint Pos, char c)
	{
	SetSS();
	const PDBChain &Chain = *m_Chain;
	const uint L = SIZE(Chain.m_Seq);
	if (Pos == 0 || Pos+1 >= L)
		return FLT_MAX;

	//vector<float> PtCA;
	//Chain.GetPt(Pos, PtCA);

	int iLo = int(Pos) - DSSParams::m_SSDensity_W;
	if (iLo < 0)
		iLo = 0;
	int iHi = int(Pos) + DSSParams::m_SSDensity_W;
	if (iHi >= int(L))
		iHi = int(L)-1;
	float D = 0;
	float Dc = 0;
	//vector<float> Pt2;
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + DSSParams::m_SSDensity_w >= Pos && Pos2 <= Pos + DSSParams::m_SSDensity_w)
			continue;
		char c2 = m_SS[Pos2];
		float Dist = Chain.GetDist(Pos, Pos2);
		float DistFactor = exp(-Dist/DSSParams::m_Density_Radius);
		D += DistFactor;
		if (c2 == c)
			Dc += DistFactor;
		}
	float r = Dc/(D + DSSParams::m_SSDensity_epsilon);
	return r;
	}

uint DSS::CalcREN(uint Pos, uint NEN) const
	{
	if (NEN == UINT_MAX)
		return UINT_MAX;

	const uint L = GetSeqLength();
	int iLo = INT_MAX;
	int iHi = INT_MAX;
	if (NEN > Pos)
		{
		iLo = int(Pos) - DSSParams::m_NEN_W;
		if (iLo < 0)
			iLo = 0;
		iHi = int(Pos) - 1;
		}
	else
		{
		iLo = int(Pos) + 1;
		iHi = int(Pos) + DSSParams::m_NEN_W;
		if (iHi >= int(L))
			iHi = int(L)-1;
		}
	if (iHi < 0)
		return UINT_MAX;
	if (iLo == INT_MAX || iHi == INT_MAX)
		return UINT_MAX;

	float MinDist = 999;
	uint MinPos = UINT_MAX;
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + DSSParams::m_NEN_w >= Pos && Pos2 <= Pos + DSSParams::m_NEN_w)
			continue;
		float Dist = m_Chain->GetDist(Pos, Pos2);
		if (Dist < MinDist)
			{
			MinDist = Dist;
			MinPos = Pos2;
			}
		}
	return MinPos;
	}

uint DSS::CalcPEN(uint Pos) const
	{
	const uint L = GetSeqLength();
	float MinDist = 999;
	uint MinPos = UINT_MAX;
	uint w = DSSParams::m_NEN_w;
	for (uint Pos2 = Pos + w; ; ++Pos2)
		{
		if (Pos2 >= L)
			break;
		float Dist = m_Chain->GetDist(Pos, Pos2);
		if (Dist < MinDist)
			{
			MinDist = Dist;
			MinPos = Pos2;
			}
		}
	return MinPos;
	}

uint DSS::CalcMEN(uint Pos) const
	{
	const uint L = GetSeqLength();
	float MinDist = 999;
	uint MinPos = UINT_MAX;
	uint w = DSSParams::m_NEN_w;
	for (int Pos2 = int(Pos) - int(w); ; --Pos2)
		{
		if (Pos2 < 0)
			break;
		float Dist = m_Chain->GetDist(Pos, Pos2);
		if (Dist < MinDist)
			{
			MinDist = Dist;
			MinPos = Pos2;
			}
		}
	return MinPos;
	}

uint DSS::CalcNEN(uint Pos) const
	{
	const uint L = GetSeqLength();
	int iLo = int(Pos) - DSSParams::m_NEN_W;
	if (iLo < 0)
		iLo = 0;
	int iHi = int(Pos) + DSSParams::m_NEN_W;
	if (iHi >= int(L))
		iHi = int(L)-1;
	float MinDist = 999;
	uint MinPos = UINT_MAX;
	for (uint Pos2 = uint(iLo); Pos2 <= uint(iHi); ++Pos2)
		{
		if (Pos2 + DSSParams::m_NEN_w >= Pos && Pos2 <= Pos + DSSParams::m_NEN_w)
			continue;
		float Dist = m_Chain->GetDist(Pos, Pos2);
		if (Dist < MinDist)
			{
			MinDist = Dist;
			MinPos = Pos2;
			}
		}
	return MinPos;
	}

uint DSS::GetNEN(uint Pos)
	{
	SetNENs();
	asserta(Pos < SIZE(m_NENs));
	return m_NENs[Pos];
	}

uint DSS::GetREN(uint Pos)
	{
	SetNENs();
	asserta(Pos < SIZE(m_RENs));
	return m_RENs[Pos];
	}

uint DSS::GetPEN(uint Pos)
	{
	SetNENs();
	asserta(Pos < SIZE(m_PENs));
	return m_PENs[Pos];
	}

uint DSS::GetMEN(uint Pos)
	{
	SetNENs();
	asserta(Pos < SIZE(m_MENs));
	return m_MENs[Pos];
	}

void DSS::SetNENs()
	{
	if (!m_NENs.empty())
		return;
	const uint L = GetSeqLength();
	asserta(SIZE(m_NENs) == 0);
	asserta(SIZE(m_RENs) == 0);
	asserta(SIZE(m_PENs) == 0);
	asserta(SIZE(m_MENs) == 0);
	m_NENs.reserve(L);
	m_RENs.reserve(L);
	m_PENs.reserve(L);
	m_MENs.reserve(L);
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		uint NEN = CalcNEN(Pos);
		uint REN = CalcREN(Pos, NEN);
		uint PEN = CalcPEN(Pos);
		uint MEN = CalcMEN(Pos);
		m_NENs.push_back(NEN);
		m_RENs.push_back(REN);
		m_PENs.push_back(PEN);
		m_MENs.push_back(MEN);
		}
	}

uint DSS::Get_NENSS(uint Pos)
	{
	SetSS();
	SetNENs();
	uint NEN = GetNEN(Pos);
	if (NEN == UINT_MAX)
		return SSCharToInt('~');
	asserta(NEN < SIZE(m_SS));
	char c = m_SS[NEN];
	return SSCharToInt(c);
	}

uint DSS::Get_RENSS(uint Pos)
	{
	SetSS();
	SetNENs();
	uint NEN = GetREN(Pos);
	if (NEN == UINT_MAX)
		return SSCharToInt('~');
	asserta(NEN < SIZE(m_SS));
	char c = m_SS[NEN];
	return SSCharToInt(c);
	}

float DSS::GetFloat_NENDist(uint Pos)
	{
	uint NEN = GetNEN(Pos);
	if (NEN == UINT_MAX)
		{
		if (opt(force_undef))
			return FLT_MAX;
		return (float) DSSParams::m_DefaultNENDist;
		}
	float d = m_Chain->GetDist(Pos, NEN);
	return d;
	}

float DSS::GetFloat_RENDist(uint Pos)
	{
	uint REN = GetREN(Pos);
	if (REN == UINT_MAX)
		return FLT_MAX;
	float d = m_Chain->GetDist(Pos, REN);
	return d;
	}

float DSS::GetFloat_PMDistDiff(uint Pos)
	{
	float pd = GetFloat_PENDist(Pos);
	float md = GetFloat_MENDist(Pos);
	if (pd == FLT_MAX || md == FLT_MAX)
		return FLT_MAX;
	float diff = pd - md;
	if (diff < -20)
		diff = -20;
	if (diff > 20)
		diff = 20;
	return diff;
	}

float DSS::GetFloat_PENDist(uint Pos)
	{
	uint PEN = GetPEN(Pos);
	if (PEN == UINT_MAX)
		return FLT_MAX;
	float d = m_Chain->GetDist(Pos, PEN);
	return d;
	}

float DSS::GetFloat_PMDist(uint Pos)
	{
	uint PEN = GetPEN(Pos);
	uint MEN = GetMEN(Pos);
	if (PEN == UINT_MAX || MEN == UINT_MAX)
		return FLT_MAX;
	float d = m_Chain->GetDist(PEN, MEN);
	return d;
	}

float DSS::GetFloat_MENDist(uint Pos)
	{
	uint MEN = GetMEN(Pos);
	if (MEN == UINT_MAX)
		return FLT_MAX;
	float d = m_Chain->GetDist(Pos, MEN);
	return d;
	}

uint DSS::Get_NormDens4(uint Pos)
	{
	uint ND = GetFeature(FEATURE_NormDens, Pos);
	if (ND == UINT_MAX)
		{
		if (opt(force_undef))
			return UINT_MAX;
		return UNDEFINED_ZERO_OVERLOAD;
		}
	asserta(ND < 16);
	return ND/4;
	}

uint DSS::Get_NENDist4(uint Pos)
	{
	uint ND = GetFeature(FEATURE_NENDist, Pos);
	if (ND == UINT_MAX)
		{
		if (opt(force_undef))
			return UINT_MAX;
		return UNDEFINED_ZERO_OVERLOAD;
		}
	asserta(ND < 16);
	return ND/4;
	}

// ADEHKNPQRST,CFILMVWY,G
uint DSS::Get_AA3(uint Pos)
	{
	const string &Seq = m_Chain->m_Seq;
	asserta(Pos < SIZE(Seq));
	char c = Seq[Pos];
	if (c == 'G')
		return 0;
	if (strchr("ADEHKNPQRST", c) != 0)
		return 1;
	if (strchr("CFILMVWY", c) != 0)
		return 2;
	return UNDEFINED_ZERO_OVERLOAD;
	}

// AHPST,CFILMVWY,DEKNQR,G
uint DSS::Get_AA4(uint Pos)
	{
	const string &Seq = m_Chain->m_Seq;
	asserta(Pos < SIZE(Seq));
	char c = Seq[Pos];
	if (c == 'G')
		return 0;
	if (strchr("AHPST", c) != 0)
		return 1;
	if (strchr("CFILMVWY", c) != 0)
		return 2;
	if (strchr("DEKNQR", c) != 0)
		return 3;
	return UNDEFINED_ZERO_OVERLOAD;
	}

void DSS::GetAaLetters(vector<byte> &Letters)
	{
	Letters.clear();
	const uint L = GetSeqLength();
	Letters.reserve(L);
	const string &Seq = m_Chain->m_Seq;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		char c = Seq[Pos];
		uint Letter = g_CharToLetterAmino[c];
		if (Letter >= 20)
			Letter = 0;
		Letters.push_back(byte(Letter));
		}
	}

void DSS::GetMuLetters(vector<byte> &Letters)
	{
	Letters.clear();
	const uint L = GetSeqLength();
	Letters.reserve(L);
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		uint Letter = Get_Mu(Pos);
		if (Letter == UINT_MAX)
			Letter = 0;
		asserta(Letter < 36);
		Letters.push_back(byte(Letter));
		}
	}

void DSS::GetProfile(vector<vector<byte> > &Profile)
	{
	Profile.clear();
	const uint L = GetSeqLength();
	const string &Seq = m_Chain->m_Seq;
	const uint FeatureCount = DSSParams::GetFeatureCount();
	Profile.reserve(FeatureCount);
	for (uint i = 0; i < FeatureCount; ++i)
		{
		vector<byte> ProfRow;
		ProfRow.reserve(L);
		FEATURE Feature = DSSParams::m_Features[i];
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			uint Letter = GetFeature(Feature, Pos);
			if (Letter == UINT_MAX)
				ProfRow.push_back(0);//@@TODO
			else
				{
				//asserta(Letter < 31);
				ProfRow.push_back(byte(Letter));
				}
			}
		Profile.push_back(ProfRow);
		}
	}

float DSS::GetFloatFeature(uint FeatureIndex, uint Pos)
	{
	switch (FeatureIndex)
		{
#define F(x)	case FEATURE_##x: return GetFloat_##x(Pos);
#include "floatfeatures.h"
#undef F
		}
	asserta(false);
	return FLT_MAX;
	}

uint DSS::GetFeature(FEATURE Feature, uint Pos)
	{
	return GetFeature(uint(Feature), Pos);
	}

uint DSS::GetFeature(uint FeatureIndex, uint Pos)
	{
	switch (FeatureIndex)
		{
		case FEATURE_AA:
		case FEATURE_B62:
			{
			char AminoChar = m_Chain->m_Seq[Pos];
			uint AminoLetter = g_CharToLetterAmino[AminoChar];
			if (AminoLetter >= 20)
				return UNDEFINED_ZERO_OVERLOAD;
			return AminoLetter;
			}

#define F(x)	case FEATURE_##x: return Get_##x(Pos);
#include "intfeatures.h"
#undef F

#define F(x)	case FEATURE_##x: \
		{ \
		float Value = GetFloat_##x(Pos); \
		uint Letter_new = DSSParams::ValueToInt_Feature(FEATURE_##x, (float) Value); \
		return Letter_new; \
		}
#include "floatfeatures.h"
#undef F

	default:
		break;
		}
	asserta(false);
	return UINT_MAX;
	}

float DSS::GetFloat_DstPrvHlx(uint Pos)
	{
	SetSSEs();
	const uint SSECount = SIZE(m_SSE_Mids);
	for (uint i = 0; i < SSECount; ++i)
		{
		if (m_SSE_cs[SSECount-i-1] != 'h')
			continue;
		uint Mid = m_SSE_Mids[i];
		if (Mid + DSSParams::m_SSE_Margin >= Pos)
			continue;
		float Dist = m_Chain->GetDist(Pos, Mid);
		return (float) Dist;
		}
	if (opt(force_undef))
		return FLT_MAX;
	return 0;
	}

float DSS::GetFloat_DstNxtHlx(uint Pos)
	{
	SetSSEs();
	const uint SSECount = SIZE(m_SSE_Mids);
	for (uint i = 0; i < SSECount; ++i)
		{
		if (m_SSE_cs[i] != 'h')
			continue;
		uint Mid = m_SSE_Mids[i];
		if (Mid <= Pos + DSSParams::m_SSE_Margin)
			continue;
		float Dist = m_Chain->GetDist(Pos, Mid);
		return Dist;
		}
	if (opt(force_undef))
		return FLT_MAX;
	return 0;
	}

uint DSS::Get_PENConf(uint Pos)
	{
	uint PEN = GetPEN(Pos);
	uint Letter = Get_Conf(PEN);
	return Letter;
	}

uint DSS::Get_MENConf(uint Pos)
	{
	uint MEN = GetMEN(Pos);
	uint Letter = Get_Conf(MEN);
	return Letter;
	}

uint DSS::Get_PENSS3(uint Pos)
	{
	SetSS();
	uint PEN = GetPEN(Pos);
	if (PEN == UINT_MAX)
		return UNDEFINED_ZERO_OVERLOAD;
	return Get_SS3(PEN);
	}

uint DSS::Get_PENSS4(uint Pos)
	{
	SetSS();
	uint PEN = GetPEN(Pos);
	if (PEN == UINT_MAX)
		return UNDEFINED_ZERO_OVERLOAD;
	return Get_SS(PEN);
	}

uint DSS::Get_MENSS3(uint Pos)
	{
	SetSS();
	uint MEN = GetMEN(Pos);
	if (MEN == UINT_MAX)
		return UNDEFINED_ZERO_OVERLOAD;
	return Get_SS3(MEN);
	}

uint DSS::Get_MENSS4(uint Pos)
	{
	SetSS();
	uint MEN = GetMEN(Pos);
	if (MEN == UINT_MAX)
		return UNDEFINED_ZERO_OVERLOAD;
	return Get_SS(MEN);
	}
