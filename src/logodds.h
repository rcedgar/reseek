#pragma once

class LogOdds
	{
public:
	uint m_AlphaSize = 0;
	vector<uint> m_BackgroundCounts;
	vector<vector<uint> > m_TrueCountMx;
	vector<float> m_Freqs;
	vector<float> m_FreqMx;

public:
	void Init(uint AlphaSize);
	void AddBackgroundLetter(uint Letter);
	void AddTruePair(uint Letter1, uint Letter2);
	void GetBackgroundFreqs(vector<float> &Freqs) const;
	void GetTrueFreqMx(vector<vector<float> > &Mx) const;
	float GetLogOddsMx(vector<vector<float> > &Mx) const;
	void GetLogOddsMxInt8(vector<vector<float> > &Mxd,
	  vector<vector<int8_t> > &Mxi, int8_t MaxAbsi) const;
	uint GetTrueTotal() const;
	void MxToSrc(FILE *f, const string &Name, 
	  const vector<vector<float> > &Mx) const;
	void MxToSrc2(FILE *f, const string &Name, 
	  const vector<vector<float> > &Mx, uint EffAlphaSize) const;
	void VecToSrc(FILE *f, const string &Name, 
	  const vector<float> &v) const;
	void GetSymbol(uint Letter, string &s) const;
	};

int8_t FloatToInt8(float x, float maxabsf, int8_t maxabsi);
