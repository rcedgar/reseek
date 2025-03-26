#pragma once

#include "dssparams.h"
#include "pdbchain.h"
#include "xdpmem.h"
#include "mx.h"
#include "userfields.h"
#include "mukmerfilter.h"
#include <mutex>

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
	vector<byte> m_MuRevA;
	void *m_ProfPara = 0;
	void *m_ProfParaRev = 0;
	MuKmerFilter m_MKF;
	float m_XDropScore = 0;
	string m_XDropPath;

	XDPMem m_Mem;
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

	float m_GlobalScore = FLT_MAX;
	string m_GlobalPath;

public:
	static mutex m_OutputLock;

public:
	static atomic<uint> m_AlnCount;
	static atomic<uint> m_SWCount;
	static atomic<uint> m_MuFilterDiscardCount;
	static atomic<uint> m_MuFilterInputCount;
	static atomic<uint> m_ParasailSaturateCount;
	static atomic<uint> m_XDropAlnCount;
	static atomic<uint> m_XDropDiscardCount;

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
	bool MuDPFilter();
	void ClearAlign();
	void AlignQueryTarget();
	void AlignQueryTarget_Global();
	void AlignQueryTarget_Trace();
	void AlignMKF();
	void PostAlignMKF();
	float GetMegaHSPScore(uint Lo_i, uint Lo_j, uint Len);
	void Align_NoAccel();
	float AlignMuQP(const vector<byte> &LettersA, const vector<byte> &LettersB);
	float AlignMuQP_Para();
	float AlignMuParaBags(const ChainBag &BagA, const ChainBag &BagB);
	float AlignMuQP_Para_Path(uint &LoA, uint &LoB, string &Path);
	float GetScorePosPair(const vector<vector<byte> > &ProfileA,
	  const vector<vector<byte> > &ProfileB, uint PosA, uint PosB) const;
	void GetPosABs(vector<uint> &PosAs, vector<uint> &PosBs) const;
	void CalcEvalue();
	void SetSMx_NoRev(const DSSParams &Params,
					  const vector<vector<byte> > &ProfileA,
					  const vector<vector<byte> > &ProfileB);
	void SetMuQP();
	void SetMuQP_Para();
	void AllocDProw(uint LB);
	const float * const *GetSMxData() const;
	float **GetSMxData();
	void AllocSMxData(uint LA, uint LB);
	void FreeSMxData();
	float GetDPScoreGivenPath(const vector<vector<byte> > &Profile1,
							  const vector<vector<byte> > &Profile2,
							  const string &Path) const;

// Up is true  if alignment is Query=A, Target=B
// Up is false if alignment is Query=B, Target=A
	void ToTsv(FILE *f, bool Up);
	void ToFasta2(FILE *f, bool Up) const;
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
	float GetNewTestStatistic(bool Top) const { return Top ? m_NewTestStatisticA : m_NewTestStatisticB; }
	float GetEvalue(bool Top) const { return Top ? m_EvalueA : m_EvalueB; }
	float GetAQ(bool Top) const { return Top ? m_QualityA : m_QualityB; }

	void GetRow(bool Up, bool Top, string &Row) const;

	float GetKabsch(float t[3], float u[3][3], bool Up) const;

	void GetRow_A(string &Row) const;
	void GetRow_B(string &Row) const;

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
	void ValidatePath() const;
	void GetCIGAR(string &CIGAR) const;

public:
	static void Stats();
	static float StaticSubstScore(void *UserData_this, uint PosA, uint PosB);
	};
