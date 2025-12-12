#pragma once

#include "myutils.h"
#include "pdbchain.h"
#include "features.h"
#include "dssparams.h"
#include "xdpmem.h"
#include "flatmx.h"

const uint UNDEFINED_ZERO_OVERLOAD = 0;

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
	vector<uint> m_PENs;
	vector<uint> m_MENs;
	string m_SS;

	vector<uint> m_SSE_Mids;
	vector<char> m_SSE_cs;
	vector<float> m_NUs;
	vector<float> m_NDs;
	vector<float> m_NXs;

	const vector<byte> *m_Seq_SSSA = 0;
	const vector<byte> *m_Seq_SSSB = 0;

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
		m_PENs.clear();
		m_MENs.clear();
		m_SSE_Mids.clear();
		m_SSE_cs.clear();
		m_NUs.clear();
		m_NDs.clear();
		m_NXs.clear();
		m_Seq_SSSA = 0;
		m_Seq_SSSB = 0;
		}

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
	uint GetPEN(uint Pos);
	uint GetMEN(uint Pos);

	uint CalcNEN(uint Pos) const;
	uint CalcREN(uint Pos, uint NEN) const;
	uint CalcPEN(uint Pos) const;
	uint CalcMEN(uint Pos) const;

	void SetNENs();
	void SetSS();
	void GetSSEs(uint MinLength, vector<uint> &Los,
	  vector<uint> &Lengths, vector<char> &cs) const;
	void SetSSEs();
	void Get_NU_ND(uint Pos, float &NU, float &ND) const;
	void Set_NU_ND_Vecs();
	void GetMuLetters(uint MuLetter, vector<uint> &Letters) const;
	float GetFloat_RENDist_ForMu(uint Pos);

public:
	static uint SSCharToInt(char c);
	static uint SSCharToInt3(char c);
	};
