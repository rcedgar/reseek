#pragma once

#include "myutils.h"
#include "pdbchain.h"
#include "features.h"
#include "dssparams.h"
#include "xdpmem.h"

const uint WILDCARD = 0;

class PDBChain;

// Discrete Structure States
class DSS
	{
public:
	const PDBChain *m_Chain = 0;
	uint m_L = 0;

	vector<double> m_Density_ScaledValues;
	vector<uint> m_NENs;
	vector<uint> m_RENs;
	string m_SS;

	int m_Density_W = 50;
	int m_Density_w = 3;
	int m_SSDensity_W = 50;
	int m_SSDensity_w = 8;
	double m_Density_Radius = 20.0;
	double m_NU_ND_Radius = 20.0;
	int m_NEN_W = 100;
	int m_NEN_w = 12;
	int m_NUDX_W = 50;
	double m_DefaultNENDist = 10.0;
	double m_SSDensity_epsilon = 1;
	uint m_SSE_MinLength = 8;
	uint m_SSE_Margin = 8;
	uint m_PMDelta = 8;

	vector<uint> m_SSE_Mids;
	vector<char> m_SSE_cs;
	vector<double> m_NUs;
	vector<double> m_NDs;
	vector<double> m_NXs;

	uint m_PatternAlphaSize1 = UINT_MAX;
	uint m_PatternAlphaSize = UINT_MAX;

private:
	const DSSParams *m_Params = 0;

public:
	void Init(const PDBChain &Chain)
		{
		m_Chain = &Chain;
		m_Density_ScaledValues.clear();
		m_SS.clear();
		m_NENs.clear();
		m_RENs.clear();
		m_SSE_Mids.clear();
		m_SSE_cs.clear();
		m_NUs.clear();
		m_NDs.clear();
		m_NXs.clear();

		if (optset_nbrw) m_NEN_w = opt_nbrw;
		if (optset_densw) m_Density_w = opt_densw;
		if (optset_ssdensw) m_SSDensity_w = opt_ssdensw;
		}

	void SetParams(const DSSParams &Params);

	uint GetSeqLength() const { return m_Chain->GetSeqLength(); }

	uint GetFeature(uint FeatureIndex, uint Pos);
	uint GetFeature(FEATURE Feature, uint Pos);
	double GetFloatFeature(uint FeatureIndex, uint Pos);

#define F(x)	uint Get_##x(uint Pos);
#include "intfeatures.h"
#undef F

#define F(x)	double GetFloat_##x(uint Pos);
#include "floatfeatures.h"
#undef F

#define F(x)	uint ValueToInt_##x(double Value) const;
#include "floatfeatures.h"
#undef F

	void GetProfile(vector<vector<byte> > &Profile);
	void GetMuLetters(vector<byte> &Letters);
	void GetMuKmers(const vector<byte> &MuLetters,
	  vector<uint> &Kmers);
	void GetMuKmerBits(const vector<uint> &Kmers,
	  vector<uint> &Bits);
	void GetMuLetters(vector<uint> &Letters);
	void SetDensity_ScaledValues();
	double GetDensity(uint Pos) const;
	double GetSSDensity(uint Pos, char c);
	uint GetNEN(uint Pos);
	uint GetREN(uint Pos);
	uint CalcNEN(uint Pos) const;
	uint CalcREN(uint Pos, uint NEN) const;

	void SetNENs();
	void SetSS();
	void GetSSEs(uint MinLength, vector<uint> &Los,
	  vector<uint> &Lengths, vector<char> &cs) const;
	void SetSSEs();
	void Get_NU_ND(uint Pos, double &NU, double &ND) const;
	void Set_NU_ND_Vecs();
	void GetMuLetters(uint MuLetter, vector<uint> &Letters) const;
	uint GetMuLetter(const vector<uint> &Letters) const;

public:
	static void GetBins(FEATURE F, vector<float> &Bins);
	static uint SSCharToInt(char c);
	static uint SSCharToInt3(char c);
	static uint ValueToInt(const vector<double> &Ts,
	  double Value);
	static uint GetAlphaSize(FEATURE F);
	};
