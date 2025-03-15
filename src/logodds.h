#pragma once

class LogOdds
	{
public:
	uint m_AlphaSize = 0;
	vector<uint> m_BackgroundCounts;
	vector<vector<uint> > m_TrueCountMx;

public:
	void Init(uint AlphaSize);
	void AddPair(uint Letter1, uint Letter2);
	void GetFreqs(vector<float> &Freqs) const;
	void GetFreqMx(vector<vector<float> > &Mx) const;
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
	float GetExpectedScore() const;
	void ToTsv(FILE *f) const;
	void FromTsv(FILE *f);
	void ValidateCounts() const;
	void ReadStringValue(FILE *f, const string &Name, string &Value);
	float ReadFloatValue(FILE *f, const string &Name, uint Idx = UINT_MAX);
	uint ReadIntValue(FILE *f, const string &Name, uint Idx = UINT_MAX);
	void ReadIntVec(FILE *f, const string &Name, uint Idx, vector<uint> &Vec);
	};

int8_t FloatToInt8(float x, float maxabsf, int8_t maxabsi);
