#pragma once

class LogOdds
	{
public:
	uint m_AlphaSize = 0;
	vector<uint> m_BackgroundCounts;
	vector<vector<uint> > m_TrueCountMx;
	vector<double> m_Freqs;
	vector<double> m_FreqMx;

public:
	void Init(uint AlphaSize);
	void AddBackgroundLetter(uint Letter);
	void AddTruePair(uint Letter1, uint Letter2);
	void GetBackgroundFreqs(vector<double> &Freqs) const;
	void GetTrueFreqMx(vector<vector<double> > &Mx) const;
	double GetLogOddsMx(vector<vector<double> > &Mx) const;
	uint GetTrueTotal() const;
	void MxToSrc(FILE *f, const string &Name, 
	  const vector<vector<double> > &Mx) const;
	void MxToSrc2(FILE *f, const string &Name, 
	  const vector<vector<double> > &Mx, uint EffAlphaSize) const;
	void VecToSrc(FILE *f, const string &Name, 
	  const vector<double> &v) const;
	void GetSymbol(uint Letter, string &s) const;
	};
