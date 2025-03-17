#include "myutils.h"
#include "logodds.h"
#include "dssparams.h"
#include "alpha.h"

int8_t FloatToInt8(float x, float maxabsf, int8_t maxabsi)
	{
	asserta(fabs(x) <= maxabsf);
	int i = int(std::round((x*maxabsi)/maxabsf));
	if (i > maxabsi)
		i = maxabsi;
	else if (i < -maxabsi)
		i = -maxabsi;
	int8_t i8 = int8_t(i);
	asserta(i8 == i);
	return i8;
	}

void LogOdds::GetSymbol(uint Letter, string &s) const
	{
	s.clear();
	if (m_AlphaSize <= 20)
		{
		s += char(g_LetterToCharAmino[Letter]);
		return;
		}
	else if (m_AlphaSize <= 26*2)
		{
		if (Letter < 26)
			{
			s += 'A' + Letter;
			return;
			}
		else if (Letter < 26*2)
			{
			s += 'a' + Letter;
			return;
			}
		else
			asserta(false);
		}
	Ps(s, "%02x", Letter);
	}

void LogOdds::Init(uint AlphaSize)
	{
	asserta(AlphaSize < 4096);
	m_AlphaSize = AlphaSize;
	m_TrueCountMx.clear();
	m_BackgroundCounts.clear();
	m_BackgroundCounts.resize(m_AlphaSize);
	m_TrueCountMx.resize(m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		m_TrueCountMx[i].resize(m_AlphaSize);
	}

void LogOdds::AddPair(uint Letter1, uint Letter2)
	{
	assert(m_AlphaSize > 0);
	if (Letter1 >= m_AlphaSize || Letter2 >= m_AlphaSize)
		return;
	m_BackgroundCounts[Letter1] += 1;
	m_BackgroundCounts[Letter2] += 1;
	m_TrueCountMx[Letter1][Letter2] += 1;
	m_TrueCountMx[Letter2][Letter1] += 1;
	}

uint LogOdds::GetTrueTotal() const
	{
	uint Total = 0;
	asserta(SIZE(m_TrueCountMx) == m_AlphaSize);
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		asserta(SIZE(m_TrueCountMx[Letter1]) == m_AlphaSize);
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			uint Count = m_TrueCountMx[Letter1][Letter2];
			Total += Count;
			}
		}
	return Total;
	}

void LogOdds::GetFreqs(vector<float> &Freqs) const
	{
	Freqs.clear();
	uint Total = 0;
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		Total += m_BackgroundCounts[Letter];

	float T = float(Total);
	float SumFreq = 0;
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		float Freq = m_BackgroundCounts[Letter]/T;
		Freqs.push_back(Freq);
		SumFreq += Freq;
		}
	asserta(feq(SumFreq, 1.0));
	}

void LogOdds::GetFreqMx(vector<vector<float> > &Mx) const
	{
	Mx.clear();
	Mx.resize(m_AlphaSize);
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		Mx[Letter].resize(m_AlphaSize);

	uint Total = GetTrueTotal();
	float T = float(Total);
	uint SumCount = 0;
	float SumFreq = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			uint Count = m_TrueCountMx[Letter1][Letter2];
			float Freq = Count/T;
			Mx[Letter1][Letter2] = Freq;
			SumCount += Count;
			SumFreq += Freq;
			}
		}
	asserta(SumCount == Total);
	asserta(feq(SumFreq, 1.0));
	}

float LogOdds::GetLogOddsMx(vector<vector<float> > &Mx) const
	{
	Mx.clear();
	Mx.resize(m_AlphaSize);
	vector<float> BackgroundFreqs;
	GetFreqs(BackgroundFreqs);
	vector<vector<float> > FreqMx;
	GetFreqMx(FreqMx);
	uint Total = GetTrueTotal();
	float SumFreq = 0;
	float ExpectedScore = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		Mx[Letter1].resize(m_AlphaSize);
		float f1 = BackgroundFreqs[Letter1];
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			float f2 = BackgroundFreqs[Letter2];
			float ObsFreq = FreqMx[Letter1][Letter2];
			float ExpectedFreq = float(f1*f2);
			if (ObsFreq == 0 || ExpectedFreq == 0)
				continue;
			float Ratio = ObsFreq/ExpectedFreq;
			float Score = log(Ratio);
			Mx[Letter1][Letter2] = Score;
			ExpectedScore += ObsFreq*Score;
			SumFreq += ObsFreq;
			}
		}
	if (!feq(SumFreq, 1.0))
		Die("LogOdds::GetLogOddsMx: SumFreq=%.6f", SumFreq);
	return ExpectedScore;
	}

float LogOdds::GetExpectedScore() const
	{
	vector<float> BackgroundFreqs;
	GetFreqs(BackgroundFreqs);
	vector<vector<float> > FreqMx;
	GetFreqMx(FreqMx);
	uint Total = GetTrueTotal();
	float SumFreq = 0;
	float ExpectedScore = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		float f1 = BackgroundFreqs[Letter1];
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			float f2 = BackgroundFreqs[Letter2];
			float ObsFreq = FreqMx[Letter1][Letter2];
			float ExpectedFreq = float(f1*f2);
			if (ObsFreq == 0 || ExpectedFreq == 0)
				continue;
			float Ratio = ObsFreq/ExpectedFreq;
			float Score = log(Ratio);
			ExpectedScore += ObsFreq*Score;
			SumFreq += ObsFreq;
			}
		}
	asserta(feq(SumFreq, 1.0));
	return ExpectedScore;
	}

void LogOdds::GetLogOddsMxInt8(vector<vector<float> > &Mxd,
  vector<vector<int8_t> > &Mxi, int8_t MaxAbsi) const
	{
	Mxi.clear();
	Mxi.resize(m_AlphaSize);
	float MaxAbs = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		for (uint Letter2 = Letter1; Letter2 < m_AlphaSize; ++Letter2)
			MaxAbs = max(MaxAbs, fabs(Mxd[Letter1][Letter2]));

	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		Mxi[Letter1].resize(m_AlphaSize);
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			float d = Mxd[Letter1][Letter2];
			int8_t i8 = FloatToInt8((float) d, (float) MaxAbs, MaxAbsi);
			Mxi[Letter1][Letter2] = i8;
			}
		}
	}

void LogOdds::VecToSrc(FILE *f, const string &Name, 
  const vector<float> &v) const
	{
	if (f == 0)
		return;
	asserta(SIZE(v) == m_AlphaSize);
	fprintf(f, "static float %s[%u] = {\n",
	  Name.c_str(), m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		fprintf(f, "\t%.4g, // %u\n", v[i], i);
	fprintf(f, "};\n");
	}

void LogOdds::MxToSrc(FILE *f, const string &Name, 
  const vector<vector<float> > &Mx) const
	{
	if (f == 0)
		return;
	asserta(SIZE(Mx) == m_AlphaSize);

	fprintf(f, "\n");
	fprintf(f, "const uint AlphaSize_%s = %u;\n", Name.c_str(), m_AlphaSize);
	fprintf(f, "static float ScoreMx_%s[%u][%u] = {\n",
	  Name.c_str(), m_AlphaSize, m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		{
		asserta(SIZE(Mx[i]) == m_AlphaSize);
		fprintf(f, "	{");
		for (uint j = 0; j < m_AlphaSize; ++j)
			{
			fprintf(f, " %7.4f", Mx[i][j]);
			if (j+1 != m_AlphaSize)
				fprintf(f, ",");
			}
		fprintf(f, "},\n");
		}
	fprintf(f, "};\n");
	}

void LogOdds::MxToSrc2(FILE *f, const string &Name, 
  const vector<vector<float> > &Mx, uint EffAlphaSize) const
	{
	if (f == 0)
		return;
	asserta(SIZE(Mx) == m_AlphaSize);
	asserta(SIZE(Mx) >= EffAlphaSize);

	fprintf(f, "SCOREMX%u_BEGIN(%s)\n",
	  EffAlphaSize, Name.c_str());
	for (uint i = 0; i < EffAlphaSize; ++i)
		{
		asserta(SIZE(Mx[i]) >= EffAlphaSize);
		fprintf(f, "SCOREMX%u_ROW(%s, %2u",
		  EffAlphaSize, Name.c_str(), i);
		for (uint j = 0; j < EffAlphaSize; ++j)
			{
			fprintf(f, ", %8.3f", Mx[i][j]);
			}
		fprintf(f, ")\n");
		}
	fprintf(f, "SCOREMX%u_END(%s)\n",
	  EffAlphaSize, Name.c_str());
	}

void LogOdds::ValidateCounts() const
	{
	asserta(SIZE(m_BackgroundCounts) == m_AlphaSize);
	asserta(SIZE(m_TrueCountMx) == m_AlphaSize);

	uint Sum1 = 0;
	for (uint i = 0; i < SIZE(m_BackgroundCounts); ++i)
		Sum1 += m_BackgroundCounts[i];

	uint Sum2 = 0;
	for (uint i = 0; i < SIZE(m_TrueCountMx); ++i)
		{
		uint RowSum = 0;
		asserta(SIZE(m_TrueCountMx[i]) == m_AlphaSize);
		for (uint j = 0; j < SIZE(m_TrueCountMx[i]); ++j)
			{
			RowSum += m_TrueCountMx[i][j];
			}
		asserta(RowSum == m_BackgroundCounts[i]);
		Sum2 += RowSum;
		}
	asserta(Sum1 == Sum2);
	}

void LogOdds::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;

	ValidateCounts();

	fprintf(f, "alpha_size\t%u\n", m_AlphaSize);
	fprintf(f, "expected_score\t%.4g\n", GetExpectedScore());
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		fprintf(f, "count\t%u\t%u",
				Letter, m_BackgroundCounts[Letter]);
		fprintf(f, "\n");
		}

	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		fprintf(f, "pair_counts\t%u", Letter);
		for(uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			fprintf(f, "\t%u", m_TrueCountMx[Letter][Letter2]);
			}
		fprintf(f, "\n");
		}
	}

void LogOdds::ReadIntVec(FILE *f, const string &Name, uint Idx, vector<uint> &Vec)
	{
	Vec.clear();
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	uint Value = UINT_MAX;
	if (Idx == UINT_MAX)
		{
		asserta(SIZE(Fields) == m_AlphaSize + 1);
		asserta(Fields[0] == Name);
		for (uint i = 0; i < m_AlphaSize; ++i)
			Vec.push_back(StrToUint(Fields[i+1]));
		}
	else
		{
		asserta(SIZE(Fields) == m_AlphaSize + 2);
		asserta(Fields[0] == Name);
		asserta(StrToUint(Fields[1]) == Idx);
		for (uint i = 0; i < m_AlphaSize; ++i)
			Vec.push_back(StrToUint(Fields[i+2]));
		}
	}

uint LogOdds::ReadIntValue(FILE *f, const string &Name, uint Idx)
	{
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	uint Value = UINT_MAX;
	if (Idx == UINT_MAX)
		{
		asserta(SIZE(Fields) == 2);
		asserta(Fields[0] == Name);
		Value = StrToUint(Fields[1]);
		}
	else
		{
		asserta(SIZE(Fields) == 3);
		asserta(Fields[0] == Name);
		asserta(StrToUint(Fields[1]) == Idx);
		Value = StrToUint(Fields[2]);
		}
	return Value;
	}

void LogOdds::ReadStringValue(FILE *f, const string &Name, string &Value)
	{
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 2);
	asserta(Fields[0] == Name);
	Value = Fields[1];
	}

float LogOdds::ReadFloatValue(FILE *f, const string &Name, uint Idx)
	{
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	float Value = FLT_MAX;
	if (Idx == UINT_MAX)
		{
		asserta(SIZE(Fields) == 2);
		asserta(Fields[0] == Name);
		Value = (float) StrToFloat(Fields[1]);
		}
	else
		{
		asserta(SIZE(Fields) == 3);
		asserta(Fields[0] == Name);
		asserta(StrToUint(Fields[1]) == Idx);
		Value = (float) StrToFloat(Fields[2]);
		}
	return Value;
	}

void LogOdds::FromTsv(FILE *f)
	{
	if (f == 0)
		return;

	m_AlphaSize = ReadIntValue(f, "alpha_size");
	float ES = ReadFloatValue(f, "expected_score");
	m_BackgroundCounts.clear();
	for (uint i = 0; i < m_AlphaSize; ++i)
		m_BackgroundCounts.push_back(ReadIntValue(f, "count", i));

	m_TrueCountMx.resize(m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		ReadIntVec(f, "pair_counts", i, m_TrueCountMx[i]);

	ValidateCounts();
	}

void LogOdds::GetExpectedScores(vector<float> &ExpectedScores) const
	{
	vector<vector<float> > ScoreMx;
	vector<float> Freqs;
	GetFreqs(Freqs);
	GetLogOddsMx(ScoreMx);

	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		float ExpectedScore = 0;
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			float Freq = Freqs[Letter];
			float Score = ScoreMx[Letter][Letter2];
			ExpectedScore += Freq*Score;
			}
		ExpectedScores.push_back(ExpectedScore);
		}
	}

uint LogOdds::GetBestDefaultLetter(uint WildcardLetter) const
	{
	vector<float> ExpectedScores;
	GetExpectedScores(ExpectedScores);
	asserta(SIZE(ExpectedScores) == m_AlphaSize);
	uint BestLetter = UINT_MAX;
	float BestExpectedScore = -FLT_MAX;
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		if (Letter == WildcardLetter)
			continue;
		float ES = ExpectedScores[Letter];
		if (ES > BestExpectedScore)
			{
			BestLetter = Letter;
			BestExpectedScore = ES;
			}
		}
	return BestLetter;
	}
