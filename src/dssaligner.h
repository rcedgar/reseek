#pragma once

#include "dssparams.h"
#include "pdbchain.h"
#include "xdpmem.h"
#include "mx.h"
#include "userfields.h"
#include "mukmerfilter.h"
#include <mutex>

enum SBSCORE
	{
	SBS_Invalid,
	SBS_Evalue,
	SBS_TS,
	SBS_Raw,
	SBS_OtherAlgoScore,	// DALI, TM-Align...
	};

#define SCORE_DIST	0

#if SCORE_DIST
#include "binner.h"
#define SCORE_BINS	100
#endif

class DSSAligner
	{
public:
	const PDBChain *m_ChainA = 0;
	const PDBChain *m_ChainB = 0;
	const vector<vector<byte> > *m_ProfileA = 0;
	const vector<vector<byte> > *m_ProfileB = 0;
	const vector<byte> *m_MuLettersA = 0;
	const vector<byte> *m_MuLettersB = 0;
	const vector<uint> *m_MuKmersA = 0;
	const vector<uint> *m_MuKmersB = 0;
	vector<const float *> m_ProfMu8;
	vector<const float *> m_ProfMu16;
	vector<const float *> m_ProfMuRev8;
	vector<const float *> m_ProfMuRev16;
	vector<byte> m_MuRevA8;
	vector<byte> m_MuRevA16;
	void *m_ProfPara8 = 0;
	void *m_ProfPara16 = 0;
	void *m_ProfParaRev8 = 0;
	void *m_ProfParaRev16 = 0;
	MuKmerFilter m_MKF;
	float m_XDropScore = 0;
	float m_XDropScoreFwd = 0;
	float m_XDropScoreBwd = 0;
	string m_XDropPath;

	XDPMem m_Mem;

	uint m_AlnDomIdx1 = UINT_MAX;
	uint m_AlnDomIdx2 = UINT_MAX;
	string m_Path;
	uint m_LoA = UINT_MAX;
	uint m_LoB = UINT_MAX;
	uint m_HiA = UINT_MAX;
	uint m_HiB = UINT_MAX;
	float m_PvalueA = FLT_MAX;
	float m_PvalueB = FLT_MAX;
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

	float m_SelfRevScoreA = FLT_MAX;
	float m_SelfRevScoreB = FLT_MAX;

	float *m_DProw = 0;
	uint m_DProwSize = 0;
	float m_AlnFwdScore = FLT_MAX;
	int m_MuFwdMinusRevScore = -999;
	int m_MuFwdScore8 = -999;
	int m_MuRevScore8 = -999;
	int m_MuFwdScore16 = -999;
	int m_MuRevScore16 = -999;

	vector<USERFIELD> m_UFs;

	float **m_SMx_Data = 0;
	float *m_SMx_Buffer = 0;
	size_t m_SMx_BufferSize = 0;
	uint m_SMx_Rows = 0;
	uint m_SMx_Cols = 0;
	bool m_DoTrace = false;//@@TODO

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

	void SetMuScore();
	bool MuFilter();
	void ClearAlign();
	void ClearAlign_ExceptMu();
	void AlignQueryTarget();
	void AlignQueryTarget_Trace();
	void Align_MuFilter(
	  const PDBChain &ChainA, const PDBChain &ChainB,
	  const vector<byte> &MuLettersA, const vector<uint> &MuKmersA,
	  const vector<byte> &MuLettersB,const vector<uint> &MuKmersB,
	  const vector<vector<byte> > &ProfileA, const vector<vector<byte> > &ProfileB);
	void AlignMKF();
	void PostAlignMKF();
	float GetMegaHSPScore(uint Lo_i, uint Lo_j, uint Len);
	void LogHSP(uint Lo_i, uint Lo_j, uint Len) const;
	void Align_NoAccel();
	void Align_QRev();

	int AlignMuQP_xx(const vector<byte> &LettersA, const vector<byte> &LettersB);
	int AlignMuQP_Para_xx();
	int AlignMuParaBags_xx(const ChainBag &BagA, const ChainBag &BagB);
	void SetMuQP_Para_xx();

	int AlignMuQP_Para8();
	int AlignMuQP_Para16();
	int AlignMuParaBags8(const ChainBag &BagA, const ChainBag &BagB);
	int AlignMuParaBags16(const ChainBag &BagA, const ChainBag &BagB);
	void SetMuQP_Para8();
	void SetMuQP_Para16();

	float GetDPScorePath(const vector<vector<byte> > &ProfileA,
	  const vector<vector<byte> > &ProfileB, uint LoA, uint LoB,
	  const string &Path) const;
	float GetScorePosPair(const vector<vector<byte> > &ProfileA,
	  const vector<vector<byte> > &ProfileB, uint PosA, uint PosB) const;
	float GetScoreSegPair(const vector<vector<byte> > &ProfileA,
	  const vector<vector<byte> > &ProfileB, uint PosA, uint PosB, uint n) const;
	uint GetU(const vector<uint> &Kmers1, const vector<uint> &Kmers2) const;
	void GetPosABs(vector<uint> &PosAs, vector<uint> &PosBs) const;
	void CalcEvalue();
	void SetSMx_QRev();
	void SetSMx_NoRev(const vector<vector<byte> > &ProfileA,
					  const vector<vector<byte> > &ProfileB);
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
	float GetSBScore(SBSCORE SBS, bool Up) const;
	//float GetAvgTestStatistic() const { return (m_TestStatisticA + m_TestStatisticB)/2; }
	float GetEvalue(bool Top) const { return Top ? m_EvalueA : m_EvalueB; }
	float GetPvalue(bool Top) const { return Top ? m_PvalueA : m_PvalueB; }
	float GetAQ(bool Top) const { return Top ? m_QualityA : m_QualityB; }
	float GetLengthCorrectedRawScore() const;

	void GetRow(bool Up, bool Top, bool Global, string &Row) const;

	float GetKabsch(double t[3], double u[3][3], bool Up) const;

	void GetRow_A(string &Row, bool Global) const;
	void GetRow_B(string &Row, bool Global) const;

	float XDropHSP(uint Loi_in, uint Loj_in, uint Len,
				   uint &Loi_out, uint &Loj_out,
				   uint &Hii_out, uint &Hij_out);

	float XDropHSP_Trace(uint Loi_in, uint Loj_in, uint Len,
				   uint &Loi_out, uint &Loj_out,
				   uint &Hii_out, uint &Hij_out);

	float GetPctId() const;
	float GetLDDT() const;
	float SubstScore(uint PosA, uint PosB);
	void AlignBags(const ChainBag &BagA,
				   const ChainBag &BagB);
	void AlignBagsMKF(const ChainBag &BagA,
				   const ChainBag &BagB);
	bool DoMKF_Bags(const ChainBag &BagA,
					const ChainBag &BagB) const;

public:
	static void Stats();
#if SCORE_DIST
	static void ReportScoreDist();
#endif
	static float StaticSubstScore(void *UserData_this, uint PosA, uint PosB);
	static float StaticSubstScore_Trace(void *UserData_this, uint PosA, uint PosB);
	};
