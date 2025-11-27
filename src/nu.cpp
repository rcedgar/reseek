#include "myutils.h"
#include "nu.h"
#include "dss.h"

extern int8_t Mu_S_ij_i8[36][36];
extern float musubstmx[36][36];
extern int8_t IntScoreMx_Mu[36][36];

void Nu::SetMu()
	{
	Clear();
	m_Features.push_back(FEATURE_SS3);
	m_Features.push_back(FEATURE_NENSS3);
	m_Features.push_back(FEATURE_RENDist4);

	m_ASs.push_back(3);
	m_ASs.push_back(3);
	m_ASs.push_back(4);

	m_Weights.push_back(1);
	m_Weights.push_back(1);
	m_Weights.push_back(1);

	m_ComponentScoreMxs.push_back(g_ScoreMxs2[FEATURE_SS3]);
	m_ComponentScoreMxs.push_back(g_ScoreMxs2[FEATURE_NENSS3]);
	m_ComponentScoreMxs.push_back(g_ScoreMxs2[FEATURE_RENDist4]);

// byte MuLetter = Letter_SS3 + Letter_NENSS3*3 + Letter_RENDist4*3*3;
	m_Axes.push_back(1);
	m_Axes.push_back(3);
	m_Axes.push_back(3*3);
	}

void Nu::SetComponents(
	const vector<FEATURE> &Fs,
	const vector<float> Weights)
	{
	const uint FeatureCount = SIZE(Fs);
	asserta(SIZE(Weights) == FeatureCount);
	Clear();
	m_Features = Fs;
	m_Weights = Weights;

	uint Axis = 1;
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = Fs[FeatureIdx];
		m_ComponentScoreMxs.push_back(DSSParams::GetScoreMx(F));
		uint AS = DSSParams::GetAlphaSize(F);
		m_ASs.push_back(AS);
		m_Axes.push_back(Axis);
		Axis *= AS;
		}
	}

void Nu::NuLetterToComponentLetters(byte NuLetter,
	vector<byte> &ComponentLetters) const
	{
	const uint FeatureCount = GetFeatureCount();
	const uint AS = GetAlphaSize();
	ComponentLetters.clear();
	ComponentLetters.resize(FeatureCount, 0);
	uint m = AS;
	for (uint k = 0; k < FeatureCount; ++k)
		{
		uint FeatureIdx = FeatureCount - k - 1;
		asserta(FeatureIdx < FeatureCount);
		byte ComponentLetter = NuLetter/m_Axes[FeatureIdx];
		ComponentLetters[FeatureIdx] = ComponentLetter;
		NuLetter -= ComponentLetter*m_Axes[FeatureIdx];
		}
	}

// byte MuLetter = Letter_SS3 + Letter_NENSS3*3 + Letter_RENDist4*3*3;
byte Nu::ComponentLettersToNuLetter(
	vector<byte> &ComponentLetters) const
	{
	uint NuLetter = 0;
	const uint FeatureCount = GetFeatureCount();
	asserta(SIZE(ComponentLetters) == FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		byte ComponentLetter = ComponentLetters[FeatureIdx];
		NuLetter += ComponentLetter*m_Axes[FeatureIdx];
		}
	byte b = byte(NuLetter);
	asserta(uint(b) == NuLetter);
	return b;
	}

uint Nu::GetAlphaSize() const
	{
	const uint FeatureCount = GetFeatureCount();
	asserta(SIZE(m_ASs) == FeatureCount);
	uint AS = 1;
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		AS *= m_ASs[FeatureIdx];
	return AS;
	}

void Nu::GetScoreMx(vector<vector<float> > &Mx) const
	{
	const uint NF = GetFeatureCount();
	const uint AS = GetAlphaSize();
	Mx.clear();
	Mx.resize(AS);
	for (uint NuLetter_i = 0; NuLetter_i < AS; ++NuLetter_i)
		{
		Mx[NuLetter_i].resize(AS);

		vector<byte> ComponentLetters_i;
		NuLetterToComponentLetters(NuLetter_i, ComponentLetters_i);

		for (uint NuLetter_j = 0; NuLetter_j < AS; ++NuLetter_j)
			{
			vector<byte> ComponentLetters_j;
			NuLetterToComponentLetters(NuLetter_j, ComponentLetters_j);
			asserta(SIZE(ComponentLetters_j) == NF);

			float Score = 0;
			for (uint FeatureIdx = 0; FeatureIdx < NF; ++FeatureIdx)
				{
				float Weight = m_Weights[FeatureIdx];
				const float * const *ComponentScoreMx = m_ComponentScoreMxs[FeatureIdx];
				uint ComponentLetter_i = ComponentLetters_i[FeatureIdx];
				uint ComponentLetter_j = ComponentLetters_j[FeatureIdx];
				Score += Weight*ComponentScoreMx[ComponentLetter_i][ComponentLetter_j];
				}
			Mx[NuLetter_i][NuLetter_j] = Score;
			}
		}
	}

void Nu::MxToSrcf(FILE *f, const string &Name, 
	const vector<vector<float> > &Mx) const
	{
	if (f == 0)
		return;
	const uint N = SIZE(Mx);
	fprintf(f, "static float %s[%u][%u] = {\n",
	  Name.c_str(), N, N);
	fprintf(f, "//	");
	for (uint i = 0; i < N; ++i)
		fprintf(f, "  %7u", i);
	fprintf(f, "\n");
	for (uint i = 0; i < N; ++i)
		{
		asserta(SIZE(Mx[i]) == N);
		fprintf(f, "	{");
		for (uint j = 0; j < N; ++j)
			{
			fprintf(f, " %7.4f", Mx[i][j]);
			if (j+1 != N)
				fprintf(f, ",");
			}
		fprintf(f, "}, // %u\n", i);
		}
	fprintf(f, "};\n");
	}

void Nu::MxToSrci(FILE *f, const string &Name, 
	const string &TypeName, uint w,
	const vector<vector<int> > &Mx) const
	{
	if (f == 0)
		return;
	const uint N = SIZE(Mx);
	fprintf(f, "static %s %s[%u][%u] = {\n",
	  TypeName.c_str(), Name.c_str(), N, N);
	fprintf(f, "//");
	for (uint i = 0; i < N; ++i)
		fprintf(f, "%*d ", w, i); // extra space for comma
	fprintf(f, "\n");
	for (uint i = 0; i < N; ++i)
		{
		asserta(SIZE(Mx[i]) == N);
		fprintf(f, " {");
		for (uint j = 0; j < N; ++j)
			{
			fprintf(f, "%*d", w, Mx[i][j]);
			if (j+1 != N)
				fprintf(f, ",");
			}
		fprintf(f, "}, // %u\n", i);
		}
	fprintf(f, "};\n");
	}

void Nu::FloatMxToIntMx(
	const vector<vector<float> > &Mxf,
	float Mul,
	vector<vector<int> > &Mxi) const
	{
	Mxi.clear();
	const uint N = SIZE(Mxf);
	Mxi.resize(N);
	for (uint i = 0; i < N; ++i)
		{
		asserta(SIZE(Mxf[i]) == N);
		for (uint j = 0; j < N; ++j)
			{
			int IntScore = int(round(Mul*Mxf[i][j]));
			Mxi[i].push_back(IntScore);
			}
		}
	}

void Nu::GetLetters(const PDBChain &Chain,
	vector<byte> &Letters)
	{
	const uint FeatureCount = GetFeatureCount();
	m_D.Init(Chain);
	const uint L = Chain.GetSeqLength();
	Letters.clear();
	vector<byte> ComponentLetters;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		ComponentLetters.clear();
		for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
			{
			FEATURE F = m_Features[FeatureIdx];
			uint ComponentLetter = m_D.GetFeature(F, Pos);
			if (ComponentLetter == UINT_MAX)
				ComponentLetters.push_back(m_ReplaceUndefWithThisLetter);
			else
				ComponentLetters.push_back(ComponentLetter);
			}
		byte Letter = ComponentLettersToNuLetter(ComponentLetters);
		Letters.push_back(Letter);
		}
	}

static void TestSetMu()
	{
	Nu A;
	A.SetMu();
	vector<byte> ComponentLetters;
	for (uint i = 0; i < 36; ++i)
		{
		A.NuLetterToComponentLetters(i, ComponentLetters);
		byte i2 = A.ComponentLettersToNuLetter(ComponentLetters);
		asserta(i2 == i);
		}
	vector<vector<float> > Mxf;
	A.GetScoreMx(Mxf);
	//A.MxToSrcf(g_fLog, "SetMu_ScoreMx", Mxf);

	for (uint i = 0; i < 36; ++i)
		for (uint j = 0; j < 36; ++j)
			asserta(feq(musubstmx[i][j], Mxf[i][j]));

	vector<vector<int> > Mxi;
	A.FloatMxToIntMx(Mxf, 1, Mxi);
	//A.MxToSrci(g_fLog, "int", "SetMu_ScoreMx1", 7, Mxi);

	for (uint i = 0; i < 36; ++i)
		for (uint j = 0; j < 36; ++j)
			asserta(IntScoreMx_Mu[i][j] == Mxi[i][j]);

	vector<vector<int> > Mxi3;
	A.FloatMxToIntMx(Mxf, 3, Mxi3);
	//A.MxToSrci(g_fLog, "int", "SetMu_ScoreMx3", 3, Mxi3);

	ProgressLog("TestSetMu() PASSED\n");
	}

static void TestSetMuComponents()
	{
	vector<FEATURE> Fs;
	vector<float> Weights;
	Fs.push_back(FEATURE_SS3);
	Fs.push_back(FEATURE_NENSS3);
	Fs.push_back(FEATURE_RENDist4);

	Weights.push_back(1);
	Weights.push_back(1);
	Weights.push_back(1);

	Nu A;
	A.SetComponents(Fs, Weights);
	vector<byte> ComponentLetters;
	for (uint i = 0; i < 36; ++i)
		{
		A.NuLetterToComponentLetters(i, ComponentLetters);
		byte i2 = A.ComponentLettersToNuLetter(ComponentLetters);
		asserta(i2 == i);
		}
	vector<vector<float> > Mxf;
	A.GetScoreMx(Mxf);
	//A.MxToSrcf(g_fLog, "SetMu_ScoreMx", Mxf);

	extern float musubstmx[36][36];
	for (uint i = 0; i < 36; ++i)
		for (uint j = 0; j < 36; ++j)
			asserta(feq(musubstmx[i][j], Mxf[i][j]));

	vector<vector<int> > Mxi;
	A.FloatMxToIntMx(Mxf, 1, Mxi);
	//A.MxToSrci(g_fLog, "int", "SetMu_ScoreMx1", 7, Mxi);

	extern int8_t IntScoreMx_Mu[36][36];
	for (uint i = 0; i < 36; ++i)
		for (uint j = 0; j < 36; ++j)
			asserta(IntScoreMx_Mu[i][j] == Mxi[i][j]);

	vector<vector<int> > Mxi3;
	A.FloatMxToIntMx(Mxf, 3, Mxi3);
	//A.MxToSrci(g_fLog, "int", "SetMu_ScoreMx3", 3, Mxi3);

	ProgressLog("TestSetMuComponents() PASSED\n");
	}

static void Match_Mu_S_ij_i8(int Idx1, int Idx2, int Idx3)
	{
	vector<FEATURE> fs;
	vector<float> Weights;
	fs.push_back(FEATURE_SS3);
	fs.push_back(FEATURE_NENSS3);
	fs.push_back(FEATURE_RENDist4);

	vector<FEATURE> Fs;
	Fs.push_back(fs[Idx1]);
	Fs.push_back(fs[Idx2]);
	Fs.push_back(fs[Idx3]);

	Weights.push_back(1);
	Weights.push_back(1);
	Weights.push_back(1);

	Nu A;
	A.SetComponents(Fs, Weights);
	vector<byte> ComponentLetters;
	for (uint i = 0; i < 36; ++i)
		{
		A.NuLetterToComponentLetters(i, ComponentLetters);
		byte i2 = A.ComponentLettersToNuLetter(ComponentLetters);
		asserta(i2 == i);
		}
	vector<vector<float> > Mxf;
	A.GetScoreMx(Mxf);

	vector<vector<int> > Mxi;
	A.FloatMxToIntMx(Mxf, 1, Mxi);
	//A.MxToSrci(g_fLog, "int", "SetMu_ScoreMx1", 7, Mxi);

	vector<vector<int> > Mxi3;
	A.FloatMxToIntMx(Mxf, 3, Mxi3);
	//A.MxToSrci(g_fLog, "int", "SetMu_ScoreMx3", 3, Mxi3);

	uint MatchCount = 0;
	for (uint i = 0; i < 36; ++i)
		for (uint j = 0; j < 36; ++j)
			if (Mxi3[i][j] == Mu_S_ij_i8[i][j])
				++MatchCount;

	ProgressLog("Match_Mu_S_ij_i8(%d, %d, %d) %u/%u matches\n",
		Idx1, Idx2, Idx3, MatchCount, 36*36);
	}

// Match_Mu_S_ij_i8(0, 1, 2) 196/1296 matches
// Match_Mu_S_ij_i8(0, 2, 1) 79/1296 matches
// Match_Mu_S_ij_i8(1, 2, 0) 71/1296 matches
// Match_Mu_S_ij_i8(1, 0, 2) 115/1296 matches
// Match_Mu_S_ij_i8(2, 0, 1) 59/1296 matches
// Match_Mu_S_ij_i8(2, 1, 0) 60/1296 matches
static void TryMatch()
	{
	Match_Mu_S_ij_i8(0, 1, 2);
	Match_Mu_S_ij_i8(0, 2, 1);
	Match_Mu_S_ij_i8(1, 0, 2);
	Match_Mu_S_ij_i8(1, 2, 0);
	Match_Mu_S_ij_i8(2, 1, 0);
	Match_Mu_S_ij_i8(2, 0, 1);
	}

static void TestChains(const string &ChainsFN)
	{
	Nu A;
	A.SetMu();

	DSS D;
	vector<PDBChain *> Chains;
	ReadChains(ChainsFN, Chains);
	const uint ChainCount = SIZE(Chains);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		vector<byte> MuLetters, NuLetters;
		D.Init(Chain);
		D.GetMuLetters(MuLetters);
		A.GetLetters(Chain, NuLetters);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			byte MuLetter = MuLetters[Pos];
			byte NuLetter = NuLetters[Pos];
			asserta(MuLetter == NuLetter);
			}
		}
	ProgressLog("TestChains() %u chains PASSED\n", ChainCount);
	}

void cmd_test_nu()
	{
	TestSetMu();
	TestSetMuComponents();
	TestChains(g_Arg1);
	}
