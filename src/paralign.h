#pragma once

#include "parasail.h"

class Paralign
	{
public:
	static parasail_matrix_t m_matrix;
	static int m_Open;
	static int m_Ext;
	static int m_SaturatedScore;
	static uint m_MaxLength;

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
	string m_Path;
	int m_LoQ = INT_MAX;
	int m_LoT = INT_MAX;

public:
	static void SetMu();
	static void SetBlosum62();
	static void SetMatrix(
		const vector<vector<int> > &ScoreMx,
		int Open, int Ext, int SaturatedScore);

public:
	void ClearResult()
		{
		m_Score = INT_MAX;
		m_Path.clear();
		m_LoQ = INT_MAX;
		m_LoT = INT_MAX;
		}

	void SetQuery(const string &LabelQ, const byte *Q, uint QL);
	void Align_ScoreOnly(const string &LabelT, const byte *T, uint TL);
	bool Align_Path(const string &LabelT, const byte *T, uint TL);
	};
