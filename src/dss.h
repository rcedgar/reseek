#pragma once

#include "myutils.h"
#include "pdbchain.h"
#include "features.h"
#include "dssparams.h"
#include "xdpmem.h"
#include "flatmx.h"
#include "undef_binning.h"
#include "timing.h"

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
	float m_SSDensity_epsilon = 1;
	uint m_SSE_MinLength = 8;
	uint m_SSE_Margin = 8;
	uint m_PMDelta = 8;

	vector<uint> m_SSE_Mids;
	vector<char> m_SSE_cs;

private:
	const DSSParams *m_Params = 0;

public:
	void Init(const PDBChain &Chain, const DSSParams &Params)
		{
		StartTimer(DSS_Init);
		m_Chain = &Chain;
		m_Density_ScaledValues.clear();
		m_SS.clear();
		m_NENs.clear();
		m_RENs.clear();
		m_SSE_Mids.clear();
		m_SSE_cs.clear();
		m_PlusNENs.clear();
		m_MinusNENs.clear();
		m_Params = &Params;
		EndTimer(DSS_Init);
		}

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
		return m_Chain->GetDist(i, j);
		}

public:
	static void Condense(const vector<float> &Values, uint AlphaSize,
						 UNDEF_BINNING UB, uint BestDefaultLetter, uint &DefaultLetter,
						 float &MinValue, float &MedValue, float &MaxValue,
						 float &UndefFreq, vector<float> &BinTs);
	static uint SSCharToInt(char c);
	static uint SSCharToInt3(char c);

	static uint ValueToInt_Feature(FEATURE F, float Value);

	static uint ValueToInt(float Value, UNDEF_BINNING UB, uint AlphaSize,
						   const vector<float> &Ts, uint DefaultLetter);

	static uint ValueToInt_Never(float Value, uint AlphaSize,
						   const vector<float> &Ts, uint DefaultLetter);

	static uint ValueToInt_OnlyZero(float Value, uint AlphaSize,
						   const vector<float> &Ts, uint DefaultLetter);

	static uint ValueToInt_ZeroOverload(float Value, uint AlphaSize,
						   const vector<float> &Ts, uint DefaultLetter);

	static uint ValueToInt_Default(float Value, uint AlphaSize,
						   const vector<float> &Ts, uint DefaultLetter);

	static uint ValueToInt_Ignore(float Value, uint AlphaSize,
						   const vector<float> &Ts, uint DefaultLetter);

	static const float *GetFreqVec(FEATURE F);
	static const float * const *GetFreqMx(FEATURE F);
	static const float * const *GetScoreMx(FEATURE F);

	static void SetFeature(FEATURE F, UNDEF_BINNING UB,
		const vector<float> &Freqs,
		const vector<vector<float> > &FreqMx,
		const vector<vector<float> > &ScoreMx,
		const vector<float> &BinTs);

	static uint GetBinThresholdCount(uint AlphaSize, UNDEF_BINNING UB);

	static uint GetAlphaSize(FEATURE F);
	static uint GetDefaultLetter(FEATURE F);
	static UNDEF_BINNING GetUB(FEATURE F);
	static const vector<float> &GetBinTs(FEATURE F);
	};
