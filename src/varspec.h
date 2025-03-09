#pragma once

class VarSpec
	{
public:
	string m_Name;
	double m_Min = DBL_MAX;
	double m_Max = DBL_MAX;
	double m_InitialDelta = DBL_MAX;
	uint m_SigFig = UINT_MAX;
	uint m_InitBinCount = UINT_MAX;	// for Latin hypercube
	double m_InitialValue = DBL_MAX;

public:
	void Init(const vector<string> &Names,
	  const vector<string> &Values);
	};
