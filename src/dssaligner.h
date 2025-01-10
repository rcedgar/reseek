#pragma once

#include "dssparams.h"
#include "pdbchain.h"
#include "xdpmem.h"
#include "mx.h"
#include "userfields.h"
#include "mukmerfilter.h"
#include <mutex>

#define SCORE_DIST	0

#if SCORE_DIST
#include "binner.h"
#define SCORE_BINS	100
#endif

class DSSAligner
	{
private:
	const DSSParams *m_Params = 0;

public:
	const PDBChain *m_ChainA = 0;
	const PDBChain *m_ChainB = 0;
	const vector<vector<byte> > *m_ProfileA = 0;
	const vector<vector<byte> > *m_ProfileB = 0;
	const vector<byte> *m_MuLettersA = 0;
	const vector<byte> *m_MuLettersB = 0;
	const vector<uint> *m_MuKmersA = 0;
	const vector<uint> *m_MuKmersB = 0;
	vector<const float *> m_ProfMu;
	vector<const float *> m_ProfMuRev;
	vector<const int8_t *> m_ProfMui;
	vector<const int8_t *> m_ProfMuRevi;
	vector<byte> m_MuRevA;
	void *m_ProfPara = 0;
	void *m_ProfParaRev = 0;
	MuKmerFilter m_MKF;
	float m_XDropScore = 0;
	string m_XDropPath;

	XDPMem m_Mem;
	//Mx<float> m_SM;

	uint m_AlnDomIdx1 = UINT_MAX;
	uint m_AlnDomIdx2 = UINT_MAX;
	string m_Path;
	uint m_LoA = UINT_MAX;
	uint m_LoB = UINT_MAX;
	uint m_HiA = UINT_MAX;
	uint m_HiB = UINT_MAX;
	float m_EvalueA = FLT_MAX;
	float m_EvalueB = FLT_MAX;
	float m_QualityA = FLT_MAX;
	float m_QualityB = FLT_MAX;
	float m_TestStatisticA = FLT_MAX;
	float m_TestStatisticB = FLT_MAX;
	float m_NewTestStatisticA = FLT_MAX;
	float m_NewTestStatisticB = FLT_MAX;

	uint m_Ids = UINT_MAX;
	uint m_Gaps = UINT_MAX;
	//float m_AlnDaliScore = FLT_MAX;
	float m_SelfRevScoreA = FLT_MAX;
	float m_SelfRevScoreB = FLT_MAX;

	float *m_DProw = 0;
	uint m_DProwSize = 0;
	float m_AlnFwdScore = FLT_MAX;

	vector<USERFIELD> m_UFs;

	float **m_SMx_Data = 0;
	float *m_SMx_Buffer = 0;
	size_t m_SMx_BufferSize = 0;
	uint m_SMx_Rows = 0;
	uint m_SMx_Cols = 0;

public:
	static mutex m_OutputLock;
	//static mutex m_StatsLock;

public:
	static atomic<uint> m_AlnCount;
	static atomic<uint> m_SWCount;
	static atomic<uint> m_MuFilterDiscardCount;
	static atomic<uint> m_MuFilterInputCount;
	static atomic<uint> m_ParasailSaturateCount;
	static atomic<uint> m_XDropAlnCount;
	static atomic<uint> m_XDropDiscardCount;
#if SCORE_DIST
	static vector<float> m_TSs;
#endif

public:
	DSSAligner();
	~DSSAligner();

public:
	void SetParams(const DSSParams &Params);

	void UnsetQuery();
	void SetQuery(
	  const PDBChain &Chain,
	  const vector<vector<byte> > *ptrProfile,
	  const vector<byte> *ptrMuLetters,
	  const vector<uint> *ptrMuKmers,
	  float SelfRevScore);
	void SetTarget(
	  const PDBChain &Chain,
	  const vector<vector<byte> > *ptrProfile,
	  const vector<byte> *ptrMuLetters,
	  const vector<uint> *ptrMuKmers,
	  float SelfRevScore);
	bool DoMKF() const;

	float GetMuScore();
	bool MuFilter();
	void ClearAlign();
	void AlignQueryTarget();
	void AlignQueryTarget_Trace();
	void Align_Test(
	  const PDBChain &ChainA, const PDBChain &ChainB,
	  const vector<byte> &MuLettersA, const vector<byte> &MuLettersB,
	  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB);
	void Align_MuFilter(
	  const PDBChain &ChainA, const PDBChain &ChainB,
	  const vector<byte> &MuLettersA, const vector<uint> &MuKmersA,
	  const vector<byte> &MuLettersB,const vector<uint> &MuKmersB,
	  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB);
	void AlignMKF();
	void PostAlignMKF();
	float GetMegaHSPScore(uint Lo_i, uint Lo_j, uint Len);
	void Align_NoAccel();
	void Align_QRev();
	float AlignMuQP(const vector<byte> &LettersA, const vector<byte> &LettersB);
	float AlignMuQP_Para();
	float AlignMuParaBags(const ChainBag &BagA, const ChainBag &BagB);
	float AlignMuQP_Para_Path(uint &LoA, uint &LoB, string &Path);
	float AlignMu_Int(const vector<byte> &LettersA, const vector<byte> &LettersB);
	float GetDPScorePath(const vector<vector<byte> > &ProfileA,
	  const vector<vector<byte> > &ProfileB, uint LoA, uint LoB,
	  const string &Path) const;
	float GetMuDPScorePath(const vector<byte> &LettersA,
	  const vector<byte> &LettersB, uint LoA, uint LoB,
	  float GapOpen, float GapExt, const string &Path) const;
	int GetMuDPScorePathInt(const vector<byte> &MuLettersA,
	  const vector<byte> &MuLettersB, uint LoA, uint LoB,
	  const string &Path) const;
	float GetScorePosPair(const vector<vector<byte> > &ProfileA,
	  const vector<vector<byte> > &ProfileB, uint PosA, uint PosB) const;
	float GetScoreSegPair(const vector<vector<byte> > &ProfileA,
	  const vector<vector<byte> > &ProfileB, uint PosA, uint PosB, uint n) const;
	uint GetU(const vector<uint> &Kmers1, const vector<uint> &Kmers2) const;
	void GetPosABs(vector<uint> &PosAs, vector<uint> &PosBs) const;
	void CalcEvalue();
	void CalcEvalue_AAOnly();
	void SetSMx_QRev();
	void SetSMx_NoRev(const DSSParams &Params,
					  const vector<vector<byte> > &ProfileA,
					  const vector<vector<byte> > &ProfileB);
	void SetMuQP();
	void SetMuQPi();
	void SetMuQP_Para();
	//void SetSMx_Mu();
	void AllocDProw(uint LB);
	const float * const *GetSMxData() const;
	float **GetSMxData();
	void AllocSMxData(uint LA, uint LB);
	void FreeSMxData();

// Up is true  if alignment is Query=A, Target=B
// Up is false if alignment is Query=B, Target=A
	void ToTsv(FILE *f, bool Up);
	void ToFasta2(FILE *f, bool Global, bool Up) const;
	void ToAln(FILE *f, bool Up) const;
	void PrettyAln(FILE *f, const PDBChain &A, const PDBChain &B,
	  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB,
	  uint LoA, uint LoB, const string &Path, float Quality, float Evalue) const;
	void WriteUserField(FILE *f, USERFIELD UF, bool Up) const;

// Top=true means fetch value for A, Top=false fetch B
	const PDBChain &GetChain(bool Top) const { return Top ? *m_ChainA : *m_ChainB; }
	const string &GetSeq(bool Top) const { return Top ? m_ChainA->m_Seq : m_ChainB->m_Seq; }
	const char *GetLabel(bool Top) const { return Top ? m_ChainA->m_Label.c_str() : m_ChainB->m_Label.c_str(); }
	uint GetQL(bool Top) const { return Top ? m_ChainA->GetSeqLength() : m_ChainB->GetSeqLength(); }
	uint GetTL(bool Top) const { return Top ? m_ChainB->GetSeqLength() : m_ChainA->GetSeqLength(); }
	uint GetLo(bool Top) const { return Top ? m_LoA : m_LoB; }
	uint GetHi(bool Top) const { return Top ? m_HiA : m_HiB; }
	uint GetL(bool Top) const { return Top ? SIZE(m_ChainA->m_Seq) : SIZE(m_ChainB->m_Seq); }
	float GetTestStatistic(bool Top) const { return Top ? m_TestStatisticA : m_TestStatisticB; }
	float GetNewTestStatistic(bool Top) const { return Top ? m_NewTestStatisticA : m_NewTestStatisticB; }
	//float GetAvgTestStatistic() const { return (m_TestStatisticA + m_TestStatisticB)/2; }
	float GetEvalue(bool Top) const { return Top ? m_EvalueA : m_EvalueB; }
	float GetAQ(bool Top) const { return Top ? m_QualityA : m_QualityB; }

	void GetRow(bool Up, bool Top, bool Global, string &Row) const;

	float GetKabsch(double t[3], double u[3][3], bool Up) const;

	void GetRow_A(string &Row, bool Global) const;
	void GetRow_B(string &Row, bool Global) const;

	float XDropHSP(uint Loi_in, uint Loj_in, uint Len,
				   uint &Loi_out, uint &Loj_out,
				   uint &Hii_out, uint &Hij_out);

	float GetPctId() const;
	float GetLDDT() const;
	float SubstScore(uint PosA, uint PosB);
	const DSSParams &GetParams() const { return *m_Params; }
	void AlignBags(const ChainBag &BagA,
				   const ChainBag &BagB);
	bool DoMKF_Bags(const ChainBag &BagA,
					const ChainBag &BagB) const;

public:
	static void Stats();
#if SCORE_DIST
	static void ReportScoreDist();
#endif
	static float StaticSubstScore(void *UserData_this, uint PosA, uint PosB);
	};
