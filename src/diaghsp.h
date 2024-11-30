#pragma once

class DiagHSP
	{
public:
	const byte *m_Q = 0;
	const byte *m_T = 0;
	int m_LQ = 0;
	int m_LT = 0;
	uint m_AS = 0;
	const short * const *m_ScoreMx = 0;
	short *m_QueryProfile = 0;

public:
	DiagHSP()
		{
		extern const short * const *Mu_S_ij_short;
		m_AS = 36;
		m_ScoreMx = Mu_S_ij_short;
		}

	void SetQ(const byte *Q, uint LQ)
		{
		m_Q = Q;
		m_LQ = LQ;
		}

	void SetT(const byte *T, uint LT)
		{
		m_T = T;
		m_LT = LT;
		}

	void SetQueryProfile();
	void FreeQueryProfile();
	int Search(int d, int &Lo, int &Len) const;
	int Search_Profile(int d, int &Lo, int &Len) const;
	int Search_Trace(int d, int &Lo, int &Len) const;
	int SearchBrute(int d, int &Lo, int &Len) const;
	int GetHSPScore(int d, int lo, int hi) const;
	int GetHSPScore_Trace(int d, int lo, int len) const;
	};
