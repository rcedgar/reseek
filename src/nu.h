#pragma once

#include "features.h"
#include "dss.h"

class Nu
	{
public:
	vector<FEATURE> m_Features;
	vector<float> m_Weights;
	vector<uint> m_ASs;
	vector<const float * const *> m_ComponentScoreMxs;
	vector<uint> m_Axes;
	DSS m_D;
	byte m_ReplaceUndefWithThisLetter = 0;

public:
	void Clear()
		{
		m_Features.clear();
		m_Weights.clear();
		m_ASs.clear();
		m_ComponentScoreMxs.clear();
		m_Axes.clear();
		}

	uint GetFeatureCount() const { return SIZE(m_Features); }
	uint GetAlphaSize() const;
	void GetScoreMx(vector<vector<float> > &Mx) const;
	void NuLetterToComponentLetters(byte NuLetter,
		vector<byte> &ComponentLetters) const;
	byte ComponentLettersToNuLetter(
		vector<byte> &ComponentLetters) const;
	void FloatMxToIntMx(
		const vector<vector<float> > &Mxf,
		float Mul,
		vector<vector<int> > &Mxi) const;
	void MxToSrcf(FILE *f, const string &Name, 
		const vector<vector<float> > &Mx) const;
	void MxToSrci(FILE *f,
		const string &TypeName,
		const string &Name, 
		uint w,
		const vector<vector<int> > &Mx) const;
	void SetMu();
	void SetComponents(
		const vector<FEATURE> &Fs,
		const vector<float> Weights);
	void GetLetters(const PDBChain &Chain,
		vector<byte> &Letters);
	};