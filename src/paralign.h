#pragma once

#include "parasail.h"
#include "xdpmem.h"

class Paralign
	{
public:
	static parasail_matrix_t m_matrix;
	static int m_Open;
	static int m_Ext;
	static int m_SaturatedScore;
	static uint m_MaxLength;
	static vector<vector<float> > m_SWFastSubstMx;

public:
	string m_LabelQ;
	string m_LabelT;
	const byte *m_Q = 0;
	const byte *m_T = 0;
	uint m_LQ = UINT_MAX;
	uint m_LT = UINT_MAX;

	parasail_profile_t *m_ProfQ = 0;
	parasail_result_t *m_result = 0;
	int m_Score = INT_MAX;
	string m_SemiGlobalPath;	// left-terminal gaps, ends in last M
	int m_LoQ = INT_MAX;		// always zero (?)
	int m_LoT = INT_MAX;		// always zero (?)

	XDPMem m_Mem;
	float m_SWFastScore = FLT_MAX;
	int m_SWFastScoreInt = INT_MAX;
	string m_SWFastPath;

public:
	static void SetMu();
	static void SetBlosum62();
	static void SetMatrix(
		const vector<vector<int> > &ScoreMx,
		int Open, int Ext, int SaturatedScore);
	static uint GetAlphaSize() { return m_matrix.size; }
	static void SetSWFastSubstMx();
	static int GetSubstScore(uint LetterQ, uint LetterT);

public:
	void ClearResult()
		{
	// Parasail
		m_Score = INT_MAX;
		m_SemiGlobalPath.clear();
		m_LoQ = INT_MAX;
		m_LoT = INT_MAX;

	// SWFast control
		m_SWFastScore = FLT_MAX;
		m_SWFastScoreInt = INT_MAX;
		m_SWFastPath.clear();
		}

	const char *GetLetterToChar() const;
	void SetQuery(const string &LabelQ, const byte *Q, uint QL);
	void Align_ScoreOnly(const string &LabelT, const byte *T, uint TL);
	bool Align_Path(const string &LabelT, const byte *T, uint TL);
	int ScoreAln(bool Trace = false) const;
	void WriteAln(FILE *f) const;
	void Align_SWFast();
	};
