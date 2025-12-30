#pragma once

#include "parasail.h"
#include "xdpmem.h"
#include "features.h"
#include <omp.h>

class Paralign
	{
public:
	static string m_SubstMxName;
	static parasail_matrix_t m_matrix;
	static int m_Open;
	static int m_Ext;
	static int m_SaturatedScore;
	static vector<vector<float> > m_SWFastSubstMx;
	static int m_Bits;

	static atomic<uint> m_Count8;
	static atomic<uint> m_Count16;
	static atomic<uint> m_Count8_rev;
	static atomic<uint> m_Count16_rev;
	static atomic<uint> m_CountSWFast;
	static atomic<uint> m_SaturatedCount;

	static bool m_GapLengthDist;
	static uint m_MaxGapLength;
	static vector<uint> m_GapLengthToCount;
	static omp_lock_t m_GapLengthLock;

public:
	string m_LabelQ;
	string m_LabelT;
	const byte *m_Q = 0;
	byte *m_Q_rev = 0;
	const byte *m_T = 0;
	uint m_LQ = UINT_MAX;
	uint m_LT = UINT_MAX;

	bool m_DoReverse = false;
	parasail_profile_t *m_ProfQ = 0;
	parasail_profile_t *m_ProfQ_rev = 0;
	parasail_result_t *m_result = 0;
	int m_Score = INT_MAX;
	int m_Score_rev = INT_MAX;
	string m_SemiGlobalPath;	// left-terminal gaps, ends in last M
	int m_LoQ = INT_MAX;		// always zero (?)
	int m_LoT = INT_MAX;		// always zero (?)

	XDPMem m_Mem;
	float m_SWFastScore = FLT_MAX;
	int m_SWFastScoreInt = INT_MAX;
	string m_SWFastPath;

public:
	static void InitGapLengthDist(uint MaxLen);
	static void LogGapLengthDist();

	static void Set_3Di();
	static void Set_Mu_S_k_i8();
	static void Set_Mu_hjmux();
	static void SetMu_musubstmx();
	static void SetMu_parasail_mu_8();
	static void SetMu_scop40_tm0_6_0_8_fa2();
	static void SetBlosum62();
	static void SetMatrix(
		const vector<vector<int> > &ScoreMx,
		int Open, int Ext, int SaturatedScore,
		bool SetSWFastMatrix = true);
	static void SetSubstMxByName(const string &Name);
	static void SetCompoundMx(
		const vector<FEATURE> &Fs, const vector<float> &Weights,
		int ScaleFactor, int Open, int Ext, int SaturatedScore);

	static uint GetAlphaSize() { return m_matrix.size; }
	static void SetSWFastSubstMx_FromParasailMx();
	static void SetSWFastSubstMx(
		const vector<vector<float> > &Mx,
		int Open, int Ext, bool DisableParasail);
	static int GetSubstScore(uint LetterQ, uint LetterT);

public:
	void ClearResult()
		{
	// Parasail
		m_Score = INT_MAX;
		m_Score_rev = INT_MAX;
		m_SemiGlobalPath.clear();
		m_LoQ = INT_MAX;
		m_LoT = INT_MAX;

	// SWFast control
		m_SWFastScore = FLT_MAX;
		m_SWFastScoreInt = INT_MAX;
		m_SWFastPath.clear();
		}

	const char *GetLetterToChar() const;
	void SetQueryProfile(const string &LabelQ, const byte *Q, uint QL);
	void SetQueryProfile_rev(const byte *Q, uint QL);
	void SetQueryNoProfile(const string &LabelQ, const byte *Q, uint QL);
	void Align_ScoreOnly(const string &LabelT, const byte *T, uint LT);
	void Align_ScoreOnly_rev(const string &LabelT, const byte *T, uint LT);
	bool Align_Path(const string &LabelT, const byte *T, uint LT);
	int ScoreAln(bool Trace = false) const;
	void WriteAln(FILE *f) const;
	void Align_SWFast(const string &LabelT, const byte *T, uint LT);
	void UpdateGapLengthDist(const string &Path);

public:
	static void LogMatrix();
	static void LogSWFastMatrix();
	static void ClearStats()
		{
		m_Count8 = 0;
		m_Count16 = 0;
		m_CountSWFast = 0;
		m_SaturatedCount = 0;
		}
	};
