#pragma once

#include "dssparams.h"
#include "pdbchain.h"
#include "xdpmem.h"
#include "dss.h"
#include "mx.h"
#include "userfields.h"
#include <mutex>

#define SCORE_DIST	0

#if SCORE_DIST
#include "binner.h"
#define SCORE_BINS	100
#endif

class DSSAligner
	{
public:
	const DSSParams *m_Params = 0;
	const PDBChain *m_ChainA = 0;
	const PDBChain *m_ChainB = 0;
	const vector<vector<byte> > *m_ProfileA = 0;
	const vector<vector<byte> > *m_ProfileB = 0;
	const vector<byte> *m_ComboLettersA = 0;
	const vector<byte> *m_ComboLettersB = 0;
	const vector<uint> *m_ComboKmerBitsA = 0;
	const vector<uint> *m_ComboKmerBitsB = 0;
	vector<const float *> m_ProfCombo;
	vector<const float *> m_ProfComboRev;
	vector<const int8_t *> m_ProfComboi;
	vector<const int8_t *> m_ProfComboRevi;
	void *m_ProfPara = 0;
	void *m_ProfParaRev = 0;

	XDPMem m_Mem;
	Mx<float> m_SMx;
	Mx<float> m_RevSMx;

	Mx<int8_t> m_SMx_Int;
	Mx<int8_t> m_RevSMx_Int;

	uint m_AlnDomIdx1 = UINT_MAX;
	uint m_AlnDomIdx2 = UINT_MAX;
	string m_PathA;
	uint m_LoA = UINT_MAX;
	uint m_LoB = UINT_MAX;
	uint m_HiA = UINT_MAX;
	uint m_HiB = UINT_MAX;
	float m_EvalueA = FLT_MAX;
	float m_EvalueB = FLT_MAX;
	float m_TestStatisticA = FLT_MAX;
	float m_TestStatisticB = FLT_MAX;

	float *m_DProw = 0;
	uint m_DProwSize = 0;
	float m_AlnFwdScore = FLT_MAX;

	float m_Query_Gumbel_mu = FLT_MAX;
	float m_Query_Gumbel_beta = FLT_MAX;

	float m_Target_Gumbel_mu = FLT_MAX;
	float m_Target_Gumbel_beta = FLT_MAX;

	vector<USERFIELD> m_UFs;

public:
	static mutex m_OutputLock;
	static mutex m_StatsLock;

public:
	static uint m_AlnCount;
	static uint m_SWCount;
	static uint m_ComboFilterCount;
	static uint m_UFilterCount;
	static uint m_ParasailSaturateCount;
#if SCORE_DIST
	static vector<float> m_TSs;
#endif

public:
	void SetQuery(
	  const PDBChain &Chain,
	  const vector<vector<byte> > *ptrProfile,
	  const vector<uint> *ptrComboKmerBits,
	  const vector<byte> *ptrComboLetters,
	  float Gumbel_mu, float Gumbel_beta);
	void SetTarget(
	  const PDBChain &Chain,
	  const vector<vector<byte> > *ptrProfile,
	  const vector<uint> *ptrComboKmerBits,
	  const vector<byte> *ptrComboLetters,
	  float Gumbel_mu, float Gumbel_beta);

	float GetComboScore();
	bool ComboFilter();
	bool UFilter();
	void AlignComboOnly();
	void AlignComboPath();
	void AlignQueryTarget();
	void Align_Test(
	  const PDBChain &ChainA, const PDBChain &ChainB,
	  const vector<byte> &ComboLettersA, const vector<byte> &ComboLettersB,
	  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB);
	void Align_ComboFilter(
	  const PDBChain &ChainA, const PDBChain &ChainB,
	  const vector<byte> &ComboLettersA, const vector<byte> &ComboLettersB,
	  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB);
	void Align_NoAccel();
	float AlignCombo(
	  const vector<byte> &LettersA, const vector<byte> &LettersB,
	  uint &LoA, uint &LoB, string &Path);
	float AlignComboQP(const vector<byte> &LettersA, const vector<byte> &LettersB);
	float AlignComboQP_Para();
	float AlignComboQP_Para_Path(uint &LoA, uint &LoB, string &Path);
	float AlignCombo_Int(const vector<byte> &LettersA, const vector<byte> &LettersB);
	float GetDPScorePath(const vector<vector<byte> > &ProfileA,
	  const vector<vector<byte> > &ProfileB, uint LoA, uint LoB,
	  const string &Path) const;
	float GetComboDPScorePath(const vector<byte> &LettersA,
	  const vector<byte> &LettersB, uint LoA, uint LoB,
	  float GapOpen, float GapExt, const string &Path) const;
	int GetComboDPScorePathInt(const vector<byte> &ComboLettersA,
	  const vector<byte> &ComboLettersB, uint LoA, uint LoB,
	  const string &Path) const;
	//float GetEvaluePath(  const PDBChain &ChainA, const PDBChain &ChainB,
	//  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB,
	//  uint LoA, uint LoB, const string &Path) const;
	float GetScorePosPair(const vector<vector<byte> > &ProfileA,
	  const vector<vector<byte> > &ProfileB, uint PosA, uint PosB) const;
	float GetScoreSegPair(const vector<vector<byte> > &ProfileA,
	  const vector<vector<byte> > &ProfileB, uint PosA, uint PosB, uint n) const;
	uint GetU(const vector<uint> &Kmers1, const vector<uint> &Kmers2) const;
	void CalcEvalue();
	void SetSMx_YesRev();
	void SetSMx_NoRev();
	void SetComboQP();
	void SetComboQPi();
	void SetComboQP_Para();
	void SetSMx_Combo();
	void SetSMx_Combo_Int();
	void AllocDProw(uint LB);
	void ToTsv(FILE *f, float MaxEvalue, bool IsA);
	void ToFasta2(FILE *f, float MaxEvalue, bool IsA);
	void ToAln(FILE *f, float MaxEvalue, bool IsA);
	float AdjustTS(float TS, float mu, float beta) const;
	void AppendUserField(string &s, USERFIELD UF, bool IsBA);
	float GetPctId() const;
	void GetRowA(string &Row) const;
	void GetRowB(string &Row) const;
	void GetRowAg(string &Row) const;
	void GetRowBg(string &Row) const;
#define x(type, name)	\
  type Get##name##A(bool IsA) { return IsA ? m_##name##A : m_##name##B; } \
  type Get##name##B(bool IsA) { return IsA ? m_##name##B : m_##name##A; }
	x(uint, Lo)
	x(uint, Hi)
	x(float, Evalue)
	x(float, TestStatistic)
#undef x

	const char *GetLabelA(bool IsA) { return IsA ? m_ChainA->m_Label.c_str() : m_ChainB->m_Label.c_str(); }
	const char *GetLabelB(bool IsA) { return IsA ? m_ChainB->m_Label.c_str() : m_ChainA->m_Label.c_str(); }

	const PDBChain &GetChainA(bool IsA) { return IsA ? *m_ChainA : *m_ChainB; }
	const PDBChain &GetChainB(bool IsA) { return IsA ? *m_ChainB : *m_ChainA; }

	const string &GetSeqA(bool IsA) { return IsA ? m_ChainA->m_Seq : m_ChainB->m_Seq; }
	const string &GetSeqB(bool IsA) { return IsA ? m_ChainB->m_Seq : m_ChainA->m_Seq; }

	uint GetLA(bool IsA) const { return IsA ? SIZE(m_ChainA->m_Seq) : SIZE(m_ChainB->m_Seq); }
	uint GetLB(bool IsA) const { return IsA ? SIZE(m_ChainB->m_Seq) : SIZE(m_ChainA->m_Seq); }

public:
	static void Stats();
#if SCORE_DIST
	static void ReportScoreDist();
#endif
	};
