#pragma once

class DiagHSP
	{
public:
	const byte *m_Q = 0;
	uint m_QL = 0;
	const byte *m_T = 0;
	uint m_TL = 0;

	uint m_QHi = 0;
	uint m_THi = 0;
	short m_Score = 0;
	uint m_Diag = UINT_MAX;

	const short * const *m_ScoreMx = 0;

public:
	void SetQ(const byte *Q, uint QL)
		{
		m_Q = Q;
		m_QL = QL;
		}

	void Search(const byte *T, uint TL, uint Diag);
	};
