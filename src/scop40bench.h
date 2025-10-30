#pragma once

#include "dbsearcher.h"
#include "dssaligner.h"
#include "xdpmem.h"
#include <map>
#include <mutex>

class SCOP40Bench : public DBSearcher
	{
public:
	SBSCORE m_SBS = SBS_Evalue;
	string m_Level = "sf";

	map<string, uint> m_LabelToChainIdx;

	map<string, uint> m_DomToChainIdx;
	vector<uint> m_DomIdxToChainIdx;

	vector<string> m_Doms;
	map<string, uint> m_DomToIdx;

	vector<string> m_SFs;
	map<string, uint> m_SFToIdx;

	vector<string> m_Folds;
	map<string, uint> m_FoldToIdx;

// Per-chain vectors [ChainIdx]
	vector<uint> m_DomIdxs;
	vector<uint> m_DomIdxToSFIdx;
	vector<uint> m_DomIdxToFoldIdx;

// Per hit vectors [HitIdx]
//   order is arbitrary (multi-threading)
	vector<float> m_Scores;
	vector<uint> m_DomIdx1s;
	vector<uint> m_DomIdx2s;

	vector<vector<uint> > m_DomIdxToHitIdxs;
	vector<uint> m_DomIdxToL;
	vector<vector<uint> > m_SFIdxToDomIdxs;
	vector<uint> m_SFSizes;
	vector<uint> m_ScoreOrder;
	vector<int> m_TFs;

	vector<float> m_ROCStepScores;
	vector<uint> m_ROCStepNTPs;
	vector<uint> m_ROCStepNFPs;
	vector<float> m_ROCStepSenss;
	vector<float> m_ROCStepEPQs;

	vector<float> m_CVESensVec;
	vector<float> m_CVEEPQVec;
	vector<float> m_CVEScoreVec;

	uint m_NT = UINT_MAX;
	uint m_NF = UINT_MAX;
	uint m_NI = UINT_MAX;
	uint m_nt_epq0_1 = UINT_MAX;
	uint m_nt_epq1 = UINT_MAX;
	uint m_nt_epq10 = UINT_MAX;

	uint m_ConsideredHitCount = UINT_MAX;
	uint m_IgnoredHitCount = UINT_MAX;

	// Area from EPQ 0.1 to 10
	float m_Area0 = FLT_MAX;

	// Area0 + (SEPQ0.1 + SEPQ1 + SEPQ10)/3
	float m_Area3 = FLT_MAX;

	FILE *m_fa2_tp = 0;
	FILE *m_fa2_fp = 0;

public:
	virtual void OnAln(DSSAligner &DA, bool Up);

public:
	void MakeSubset(SCOP40Bench &Subset, uint Pct) const;
	void ReadLookup(const string &FileName);
	void ClearHitsAndResults();
	float GetVeryBadScore() const;
	float GetVeryGoodScore() const;
	bool ScoreIsBetter(float Score1, float Score2) const;
	bool HitIsBetter(uint HitIdx1, uint HitIdx2) const;
	uint GetHitCount() const { return SIZE(m_Scores); }
	uint GetDomCount() const { return SIZE(m_DomIdxToSFIdx); }
	uint GetSFCount() const { return SIZE(m_SFs); }
	uint GetFoldCount() const { return SIZE(m_Folds); }
	void AddDom(const string &Dom, const string &Fold, const string &SF,
	  uint ChainIndex = UINT_MAX);
	void BuildDomSFIndexesFromDBChainLabels();
	uint GetDomIdx(const string &Dom_or_DomSlashId, bool FailOnErr = true) const;
	const PDBChain &GetChainByDomIdx(uint DomIdx) const;
	const vector<vector<byte> > &GetProfileByDomIdx(uint DomIdx) const;

	void SetSFIdxToDomIdxs();
	void SetDomIdxToL();
	void ReadHits(const string &FN);
	void WriteBit(const string &FileName) const;
	void ReadBit(const string &FileName);
	int IsT(uint DomIdx1, uint DomIdx2) const;
	void LoadHitsFromTsv(const string &FileName);
	float GetMeanLength(uint SFIdx) const;
	void StoreScore(uint ChainIdx1, uint ChainIdx2, float Score12);

	void WriteFasta2s(DSSAligner &DA) const;

// ROC analysis
	void SetStats(float MaxFPR);
	void SetTFs();
	void SetNXs();
	void SetScoreOrder();
	//void SetTSOrder();
	
	void SetROCSteps();
	void SetCVE();

	void WriteCVE(FILE *f) const;

	uint GetNTPAtEPQThreshold(const vector<uint> &NTPs,
	  const vector<uint> &NFPs, float EPQ) const;
	float GetTPRAtEPQThreshold(const vector<uint> &NTPs,
	  const vector<uint> &NFPs, float EPQ) const;
	float GetEPQAtEvalueThreshold(const vector<float> &Evalues,
	  const vector<uint> &NFPs, float Evalue) const;
	float GetEvalueAtEPQThreshold(const vector<float> &Evalues,
	  const vector<uint> &NFPs, float NFP) const;

	float AlignDomPair(uint ThreadIndex, uint Dom1, uint Dom2,
	  uint &Lo1, uint &Lo2, string &Path);
	void WriteSummary();
	void WriteOutput();
	void LogFirstFewDoms() const;
	void LogFirstFewHits() const;
	void WriteSteps(const string &FN, bool WithHdr = false) const;
	void WriteSortedHits(const string &FN) const;
	void RoundScores();
	void SetArea();

public:
	virtual bool Reject(DSSAligner &DA, bool Up) const { return false; }
	virtual void OnSetup();

public:
	static void ParseScopLabel(const string &Label, string &Dom,
	  string &Cls, string &Fold, string &SF, string &Fmy,
	  bool MissingOk = false);
	static void GetDomFromLabel(const string &Label, string &Dom);
	static void GetDomSFFromLabel(const string &Label, string &Dom, string &SF);
	static bool IsTP_SF(const string &Label1, const string &Label2);
	static const char *SBSToStr(SBSCORE SBS)
		{
		switch (SBS)
			{
			case SBS_Evalue:	return "Evalue";
			case SBS_TS:		return "TS";
			case SBS_Raw:		return "Raw";
			case SBS_OtherAlgoScore:	return "OtherAlgo";
			}
		Die("SBSToStr(%d)", int(SBS));
		return "*ERROR(";
		}
	};
