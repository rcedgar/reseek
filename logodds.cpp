#include "myutils.h"
#include "logodds.h"
#include "dssparams.h"
#include "alpha.h"

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
	m_AlphaSize = AlphaSize;
	m_TrueCountMx.clear();
	m_Freqs.clear();
	m_FreqMx.clear();
	m_BackgroundCounts.resize(m_AlphaSize);
	m_TrueCountMx.resize(m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		m_TrueCountMx[i].resize(m_AlphaSize);
	}

void LogOdds::AddBackgroundLetter(uint Letter)
	{
	assert(m_AlphaSize > 0);
	if (Letter >= m_AlphaSize)
		return;
	asserta(SIZE(m_BackgroundCounts) == m_AlphaSize);
	m_BackgroundCounts[Letter] += 1;
	}

void LogOdds::AddTruePair(uint Letter1, uint Letter2)
	{
	assert(m_AlphaSize > 0);
	if (Letter1 >= m_AlphaSize || Letter2 >= m_AlphaSize)
		return;
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

void LogOdds::GetBackgroundFreqs(vector<double> &Freqs) const
	{
	Freqs.clear();
	uint Total = 0;
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		Total += m_BackgroundCounts[Letter];

	double T = float(Total);
	double SumFreq = 0;
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		double Freq = m_BackgroundCounts[Letter]/T;
		Freqs.push_back(Freq);
		SumFreq += Freq;
		}
	asserta(feq(SumFreq, 1.0));
	}

void LogOdds::GetTrueFreqMx(vector<vector<double> > &Mx) const
	{
	Mx.clear();
	Mx.resize(m_AlphaSize);
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		Mx[Letter].resize(m_AlphaSize);

	uint Total = GetTrueTotal();
	double T = double(Total);
	uint SumCount = 0;
	double SumFreq = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			uint Count = m_TrueCountMx[Letter1][Letter2];
			double Freq = Count/T;
			Mx[Letter1][Letter2] = Freq;
			SumCount += Count;
			SumFreq += Freq;
			}
		}
	asserta(SumCount == Total);
	asserta(feq(SumFreq, 1.0));
	}

double LogOdds::GetLogOddsMx(vector<vector<double> > &Mx) const
	{
	Mx.clear();
	Mx.resize(m_AlphaSize);
	vector<double> BackgroundFreqs;
	GetBackgroundFreqs(BackgroundFreqs);
	vector<vector<double> > FreqMx;
	GetTrueFreqMx(FreqMx);
	uint Total = GetTrueTotal();
	double SumFreq = 0;
	double ExpectedScore = 0;
	for (uint Letter1 = 0; Letter1 < m_AlphaSize; ++Letter1)
		{
		Mx[Letter1].resize(m_AlphaSize);
		double f1 = BackgroundFreqs[Letter1];
		for (uint Letter2 = 0; Letter2 < m_AlphaSize; ++Letter2)
			{
			double f2 = BackgroundFreqs[Letter2];
			double ObsFreq = FreqMx[Letter1][Letter2];
			double ExpectedFreq = double(f1*f2);
			if (ObsFreq == 0 || ExpectedFreq == 0)
				continue;
			double Ratio = ObsFreq/ExpectedFreq;
			double Score = log(Ratio);
			Mx[Letter1][Letter2] = Score;
			ExpectedScore += ObsFreq*Score;
			SumFreq += ObsFreq;
			}
		}
	asserta(feq(SumFreq, 1.0));
	return ExpectedScore;
	}

void LogOdds::VecToSrc(FILE *f, const string &Name, 
  const vector<double> &v) const
	{
	if (f == 0)
		return;
	asserta(SIZE(v) == m_AlphaSize);
	fprintf(f, "static double %s[%u] = {\n",
	  Name.c_str(), m_AlphaSize);
	for (uint i = 0; i < m_AlphaSize; ++i)
		fprintf(f, "\t%.4g, // %u\n", v[i], i);
	fprintf(f, "};\n");
	}

void LogOdds::MxToSrc(FILE *f, const string &Name, 
  const vector<vector<double> > &Mx) const
	{
	if (f == 0)
		return;
	asserta(SIZE(Mx) == m_AlphaSize);

	fprintf(f, "\n");
	fprintf(f, "const uint AlphaSize_%s = %u;\n", Name.c_str(), m_AlphaSize);
	fprintf(f, "static double ScoreMx_%s[%u][%u] = {\n",
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
  const vector<vector<double> > &Mx, uint EffAlphaSize) const
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
