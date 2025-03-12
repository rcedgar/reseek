#pragma once

#include "myutils.h"
#include "pdbchain.h"
#include "features.h"
#include "dssparams.h"
#include "xdpmem.h"
#include "flatmx.h"

class PDBChain;
class DSSAligner;

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
	int m_NEN_W = 100;
	int m_NEN_w = 12;
	int m_NUDX_W = 50;
	float m_SSDensity_epsilon = 1;
	uint m_SSE_MinLength = 8;
	uint m_SSE_Margin = 8;
	uint m_PMDelta = 8;

	vector<uint> m_SSE_Mids;
	vector<char> m_SSE_cs;
#if CACHE_DIST_MAX
	const FlatMx<float> *m_DistMx = 0;
#endif
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
		m_PlusNENs.clear();
		m_MinusNENs.clear();

#if CACHE_DIST_MAX
		asserta(Chain.m_DistMx != 0);
		m_DistMx = Chain.m_DistMx;
#endif
		}

	void SetParams(const DSSParams &Params);

	uint GetSeqLength() const { return m_Chain->GetSeqLength(); }

	uint GetFeature(uint FeatureIndex, uint Pos);
	uint GetFeature(FEATURE Feature, uint Pos);
	uint GetFeatureLo(FEATURE Feature, uint Pos);

// Return FLT_MAX if undefined
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
	uint GetPlusNEN(uint Pos);
	uint GetMinusNEN(uint Pos);

	uint CalcREN(uint Pos, uint NEN) const;
	uint CalcNEN(uint Pos) const;
	uint CalcPlusNEN(uint Pos) const;
	uint CalcMinusNEN(uint Pos) const;

	void SetNENs();
	void SetSS();
	void GetSSEs(uint MinLength, vector<uint> &Los,
	  vector<uint> &Lengths, vector<char> &cs) const;
	void SetSSEs();
	void GetMuLetters(uint MuLetter, vector<uint> &Letters) const;
	uint GetMuLetter(const vector<uint> &Letters) const;

	float GetDist(uint i, uint j) const
		{
#if CACHE_DIST_MAX
		assert(m_DistMx != 0);
		return m_DistMx->Get(i, j);
#else
		return m_Chain->GetDist(i, j);
#endif
		}

public:
	static void Condense(const vector<float> &Values, uint AlphaSize,
						 float &MinValue, float &MedValue,
						 float &MaxValue, float &UndefFreq,
						 vector<float> &BinTs);
	static uint SSCharToInt(char c);
	static uint SSCharToInt3(char c);
	static uint ValueToInt(const vector<float> &Ts, float Value);
	static uint GetAlphaSize(FEATURE F);
	};

float GetSelfRevScore(DSSAligner &DA,
	const DSSParams &Params,
	const PDBChain &Chain,
	const vector<vector<byte> > &Profile,
	const vector<byte> *ptrMuLetters,
	const vector<uint> *ptrMuKmers);

const float *GetFreqVec(FEATURE F);
const float * const *GetFreqMx(FEATURE F);
const float * const *GetScoreMx(FEATURE F);
uint GetAlphaSize(FEATURE F);
