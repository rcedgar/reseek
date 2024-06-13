#pragma once

#include "dbsearcher.h"
#include "dssaligner.h"
#include "xdpmem.h"
#include <map>
#include <mutex>

class SCOP40Bench : public DBSearcher
	{
public:
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

	bool m_ScoresAreEvalues = opt_scores_are_evalues;

// Per hit vectors [HitIdx]
//   order is arbitrary (multi-threading)
	vector<float> m_Scores;
	vector<uint> m_DomIdx1s;
	vector<uint> m_DomIdx2s;

	uint m_ChainIndex1 = UINT_MAX;
	uint m_ChainIndex2 = UINT_MAX;

	vector<float> m_DomIdxToScoreFirstFP;
	vector<uint> m_DomIdxToHitIdxFirstFP;
	vector<vector<uint> > m_DomIdxToHitIdxs;
	vector<uint> m_DomIdxToL;
	vector<uint> m_DomIdxToSens1FP;
	vector<uint> m_DomIdxToTP1Count;
	vector<vector<uint> > m_SFIdxToDomIdxs;
	vector<uint> m_SFSizes;
	vector<uint> m_ScoreOrder;
	vector<int> m_TFs;

	vector<float> m_ROCStepScores;
	vector<uint> m_ROCStepNTPs;
	vector<uint> m_ROCStepNFPs;
	vector<float> m_SmoothTPRs;
	vector<float> m_SmoothFPRs;
	vector<float> m_SmoothScores;
	vector<uint> m_SmoothNTPs;
	vector<uint> m_SmoothNFPs;
	uint m_NT = UINT_MAX;
	uint m_NF = UINT_MAX;
	uint m_NI = UINT_MAX;
	uint m_nt_epq0_1 = UINT_MAX;
	uint m_nt_epq1 = UINT_MAX;
	uint m_nt_epq10 = UINT_MAX;
	uint m_nt_firstfp = UINT_MAX;

	string m_Level;
	uint m_ConsideredHitCount = UINT_MAX;
	uint m_IgnoredHitCount = UINT_MAX;

public:
	virtual void OnSetup();
	virtual void OnAln(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA);
	virtual void OnAlnBA(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA);

public:
	void ReadLookup(const string &FileName);
	void ClearHits();
	float GetVeryBadScore() const;
	bool KeepScore(float Score) const;
	bool ScoreIsBetter(float Score1, float Score2) const;
	bool HitIsBetter(uint HitIdx1, uint HitIdx2) const;
	uint GetHitCount() const { return SIZE(m_Scores); }
	uint GetDomCount() const { return SIZE(m_DomIdxToSFIdx); }
	uint GetSFCount() const { return SIZE(m_SFs); }
	uint GetFoldCount() const { return SIZE(m_Folds); }
	void AddDom(const string &Dom, const string &Fold, const string &SF,
	  uint ChainIndex = UINT_MAX);
	void BuildDomSFIndexesFromQueryChainLabels();
	uint GetDomIdx(const string &Dom_or_DomSlashId, bool FailOnErr = true) const;
	const PDBChain &GetChainByDomIdx(uint DomIdx) const;
	const vector<vector<byte> > &GetProfileByDomIdx(uint DomIdx) const;

	void ScanDomHits();
	void SetDomIdxToHitIdxs();
	void SetSFIdxToDomIdxs();
	void SetDomIdxToL();
	uint GetSens1stFP();
	void GetTPs1FP(vector<uint> &Doms1, vector<uint> &Doms2);
	void ReadHits(const string &FN);
	void WriteBit(const string &FileName) const;
	void ReadBit(const string &FileName);
	int IsT(uint DomIdx1, uint DomIdx2) const;
	void LoadHitsFromTsv(const string &FileName);
	float GetMeanLength(uint SFIdx) const;
	void StoreScore(uint ChainIdx, uint ChainIdx2, float Score12);

// ROC analysis
	void SetStats(float MaxFPR);
	void SetTFs();
	void SetNXs();
	void SetScoreOrder();
	bool SmoothROCSteps(const vector<float> &Scores,
	  const vector<uint> &NTPs, const vector<uint> &NFPs,
	  uint N, float MaxFPR, vector<float> &SScores,
	  vector<uint> &SNTPs, vector<uint> &SNPS,
	  vector<float> &STPRs, vector<float> &SFPRs) const;
	void GetROCSteps(vector<float> &ScoreSteps,
	   vector<uint> &NTPs, vector<uint> &NFPs);
	void WriteSensVsErr(FILE *f, uint N);
	void ROCToTsv(const string &FileName, float MaxFPR);
	uint GetNTPAtEPQThreshold(const vector<uint> &NTPs,
	  const vector<uint> &NFPs, float EPQ) const;
	float GetTPRAtEPQThreshold(const vector<uint> &NTPs,
	  const vector<uint> &NFPs, float EPQ) const;
	float GetEPQAtEvalueThreshold(const vector<float> &Evalues,
	  const vector<uint> &NFPs, float Evalue) const;
	float GetEvalueAtEPQThreshold(const vector<float> &Evalues,
	  const vector<uint> &NFPs, float NFP) const;
	void ROCStepsToTsv(const string &FileName,
	   vector<float> &ScoreSteps,
	   vector<uint> &NTPs, vector<uint> &NFPs,
	   vector<float> &TPRs, vector<float> &FPRs) const;
	float AlignDomPair(uint ThreadIndex, uint Dom1, uint Dom2,
	  uint &Lo1, uint &Lo2, string &Path);
	void WriteSummary();
	void WriteOutput();
	void LogSens1FPReport();

public:
	static void ParseScopLabel(const string &Label, string &Dom,
	  string &Cls, string &Fold, string &SF, string &Fmy,
	  bool MissingOk = false);
	static void GetDomFromLabel(const string &Label, string &Dom);
	};
