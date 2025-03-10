#pragma once

#include "myutils.h"
#include "pdbchain.h"
#include "features.h"
#include "dssparams.h"
#include "xdpmem.h"
#include "flatmx.h"

const uint WILDCARD = 0;

class PDBChain;

// Discrete Structure States
class DSS
	{
public:
	const PDBChain *m_Chain = 0;
	uint m_L = 0;

	vector<float> m_Density_ScaledValues;
	vector<uint> m_NENs;
	vector<uint> m_RENs;
	vector<uint> m_PlusNENs;
	vector<uint> m_MinusNENs;
	string m_SS;

	int m_Density_W = 50;
	int m_Density_w = 3;
	int m_SSDensity_W = 50;
	int m_SSDensity_w = 8;
	float m_Density_Radius = 20.0;
	float m_NU_ND_Radius = 20.0;
	int m_NEN_W = 100;
	int m_NEN_w = 12;
	int m_NUDX_W = 50;
	float m_DefaultNENDist = 10.0;
	float m_SSDensity_epsilon = 1;
	uint m_SSE_MinLength = 8;
	uint m_SSE_Margin = 8;
	uint m_PMDelta = 8;

	vector<uint> m_SSE_Mids;
	vector<char> m_SSE_cs;
	vector<float> m_NUs;
	vector<float> m_NDs;
	vector<float> m_NXs;

	const FlatMx<float> *m_DistMx = 0;

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

		if (optset_nbrw) m_NEN_w = opt(nbrw);
		if (optset_densw) m_Density_w = opt(densw);
		if (optset_ssdensw) m_SSDensity_w = opt(ssdensw);

		asserta(Chain.m_DistMx != 0);
		m_DistMx = Chain.m_DistMx;
		}

	void SetParams(const DSSParams &Params);

	uint GetSeqLength() const { return m_Chain->GetSeqLength(); }

	uint GetFeature(uint FeatureIndex, uint Pos);
	uint GetFeature(FEATURE Feature, uint Pos);
	float GetFloatFeature(uint FeatureIndex, uint Pos);

#define F(x)	uint Get_##x(uint Pos);
#include "intfeatures.h"
#undef F

#define F(x)	float GetFloat_##x(uint Pos);
#include "floatfeatures.h"
#undef F

#define F(x)	uint ValueToInt_##x(float Value) const;
#include "floatfeatures.h"
#undef F

	void GetProfile(vector<vector<byte> > &Profile);
	void GetMuLetters(vector<byte> &Letters);
	void GetAaLetters(vector<byte> &Letters);
	void GetMuKmers(const vector<byte> &MuLetters,
	  vector<uint> &Kmers, const string &PatternStr);
	void GetMuLetters(vector<uint> &Letters);
	void SetDensity_ScaledValues();
	float GetDensity(uint Pos) const;
	float GetSSDensity(uint Pos, char c);
	uint GetNEN(uint Pos);
	uint GetREN(uint Pos);
	uint CalcNEN(uint Pos) const;
	uint CalcREN(uint Pos, uint NEN) const;
	void SetNENs();
	void SetSS();
	void GetSSEs(uint MinLength, vector<uint> &Los,
	  vector<uint> &Lengths, vector<char> &cs) const;
	void SetSSEs();
	void Get_NU_ND(uint Pos, float &NU, float &ND) const;
	void Set_NU_ND_Vecs();
	void GetMuLetters(uint MuLetter, vector<uint> &Letters) const;
	uint GetMuLetter(const vector<uint> &Letters) const;

	float GetDist(uint i, uint j) const
		{
		assert(m_DistMx != 0);
		return m_DistMx->Get(i, j);
		}

public:
	static uint SSCharToInt(char c);
	static uint SSCharToInt3(char c);
	static uint ValueToInt(const vector<float> &Ts,
	  float Value);
	static uint GetAlphaSize(FEATURE F);
	};
