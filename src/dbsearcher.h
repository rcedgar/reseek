#pragma once

#include "profileloader.h"
#include "dbsearcher.h"
#include "dssaligner.h"
#include "xdpmem.h"
#include <atomic>
#include <map>
#include <mutex>
#include "mukmerfilter.h"
#include "triangle.h"

class ChainReader2;

class DBSearcher
	{
public:
	~DBSearcher();

public:
	mutex m_Lock;
	uint m_ThreadCount = UINT_MAX;
	vector<DSSAligner *> m_DAs;
	vector<XDPMem *> m_Mems;

	vector<PDBChain *> m_DBChains;
	bool m_QuerySelf = false;

// Per-chain vectors [ChainIdx]
	vector<vector<vector<byte> > *> m_DBProfiles;
	vector<vector<byte> *> m_DBMuLettersVec;
	vector<vector<uint> *> m_DBMuKmersVec;
	vector<float> m_DBSelfRevScores;

	atomic<uint> m_PairIndex = UINT_MAX;
	uint m_PairCount = UINT_MAX;
	uint m_NextChainIndex1 = UINT_MAX;
	uint m_NextChainIndex2 = UINT_MAX;
	uint m_NextQueryIdx = UINT_MAX;
	uint m_NextDBIdx = UINT_MAX;

	atomic<uint> m_ProcessedQueryCount = 0;
	atomic<uint> m_ProcessedPairCount = 0;
	atomic<uint> m_HitCount = 0;
	atomic<uint> m_QPCacheHits = 0;
	atomic<uint> m_QPCacheMisses = 0;

	uint m_FilterRejects = 0;
	uint m_XAlignCount = 0;
	uint m_SWAlignCount = 0;
	uint m_UFilterCount = 0;
	double m_MinEvalue = -1;
	double m_MaxEvalue = 10;
	uint m_Secs = UINT_MAX;
	float m_AlnsPerThreadPerSec = FLT_MAX;
	time_t m_LastProgress = 0;

	bool m_RecalcSelfRevScores = false;

	uint8_t *m_dope = 0;
	uint32_t m_dope_nhit = 0;
	uint32_t *m_dope_ks = 0;


public:
	void Setup();
	void InitEmpty();
	void ClearStats();
	void LoadDB(const string &DBFN);
	uint GetDBChainCount() const { return SIZE(m_DBChains); }
	void SetSelfRevScores();

	void RunQuery(ChainReader2 &QCR);
	void RunSelf(bool ShowStats = true);

	void ThreadBodyQuery(uint ThreadIndex, ChainReader2 *ptrQueryCR);
	void ThreadBodySelf(uint ThreadIndex);
	void ThreadBodySelf_MuFilterOnly(uint ThreadIndex);

	uint GetDBSize() const;
	bool GetNextPairSelf(uint &ChainIndex1, uint &ChainIndex2);
	void RunStats() const;
	void AddChain(PDBChain *ptrChain, vector<vector<byte> > *ptrProfile,
	  vector<byte> *ptrMuLetters);
	void ShuffleProfiles();
	void ShuffleProfile(vector<vector<byte> > &Profile);

	void ReadDope(const string &FN);
	bool in_dope(uint k) const
		{
		if (m_dope == 0) return true;
		byte b = m_dope[k/8];
		return b & (1 << k%8);
		}

	bool in_dope(uint i, uint j) const
		{
		if (m_dope == 0) return true;
		uint k = triangle_ij_to_k2(i, j, SIZE(m_DBChains));
		return in_dope(k);
		}

	virtual bool in_dope_labels(const string &label_i, const string &label_j) const
		{
		return false;
		}

	float get_missing_score() const
		{
		if (opt(scores_are_evalues))
			return 9999;
		else
			return -9999;
		}

public:
	virtual bool Reject(DSSAligner &DA, bool Up) const;

public:
	virtual void OnSetup() {}
	void BaseOnAln(DSSAligner &DA, bool Up);
	virtual void OnAln(DSSAligner &DA, bool Up) {}

public:
	static void StaticThreadBodyQuery(uint ThreadIndex, DBSearcher *ptrDBS, ChainReader2 *ptrQueryCR);
	static void StaticThreadBodySelf(uint ThreadIndex, DBSearcher *ptrDBS);
	};
