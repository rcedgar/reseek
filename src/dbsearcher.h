#pragma once

#include "dbsearcher.h"
#include "dssaligner.h"
#include "xdpmem.h"
#include <atomic>
#include <map>
#include <mutex>

class DBSearcher
	{
public:
	const DSSParams *m_Params = 0;
	uint m_ThreadCount = UINT_MAX;
	vector<DSSAligner *> m_DAs;
	vector<XDPMem *> m_Mems;
	DSS m_D;

// m_Chains has m_DBChainCount + m_QueryChainCount = 
//   m_ChainCount chain pointers
	vector<PDBChain *> m_Chains;
	vector<PDBChain *> m_QueryChains;
	vector<PDBChain *> m_DBChains;
	uint m_DBChainCount = UINT_MAX;
	uint m_QueryChainCount = UINT_MAX;
	bool m_QuerySelf = false;

// Per-chain vectors [ChainIdx]
	vector<vector<vector<byte> > > m_Profiles;
	vector<vector<byte> > m_ComboLettersVec;
	vector<vector<uint> > m_KmersVec;
	vector<vector<uint> > m_KmerBitsVec;

	mutex m_Lock;
	uint m_PairIndex = UINT_MAX;
	uint m_PairCount = UINT_MAX;
	uint m_NextChainIndex1 = UINT_MAX;
	uint m_NextChainIndex2 = UINT_MAX;
	uint m_NextQueryIdx = UINT_MAX;
	uint m_NextDBIdx = UINT_MAX;
	uint m_ChainCount = UINT_MAX;
	uint m_ProcessedQueryCount = UINT_MAX;
	atomic<uint> m_QPCacheHits;
	atomic<uint> m_QPCacheMisses;

	uint m_FilterRejects = 0;
	uint m_XAlignCount = 0;
	uint m_SWAlignCount = 0;
	uint m_UFilterCount = 0;
	float m_MaxEvalue = FLT_MAX;
	FILE *m_fTsv = 0;
	FILE *m_fAln = 0;
	FILE *m_fFasta2 = 0;
	uint m_Secs = UINT_MAX;
	float m_AlnsPerThreadPerSec = FLT_MAX;

// Calibration
	bool m_CollectTestStats = false;
	vector<vector<float> > m_TestStatsVec;

public:
	const PDBChain &GetQueryChain(uint Idx) const;
	const char *GetQueryLabel(uint Idx) const;
	void Setup(const DSSParams &Params);
	uint GetProfileCount() const { return SIZE(m_Profiles); }
	uint GetChainCount() const { return SIZE(m_Chains); }
	bool GetNextPair(uint &ChainIndex1, uint &ChainIndex2);
	bool GetNextPairQuerySelf(uint &ChainIndex1, uint &ChainIndex2);
	const PDBChain &GetChain(uint ChainIndex) const;
	void ReadChains(const string &QueryCalFileName, 
	  const string &DBCalFileName = "");
	void SetProfiles();
	void SetKmersVec();
	void Run();
	void RunUSort();
	uint GetDBSize() const;
	void USort(uint QueryChainIndex, vector<uint> &Idxs,
	  vector<uint> &DBChainIndexes);
	uint GetDBChainIndex(uint Idx) const;
	void Thread(uint ThreadIndex);
	void ThreadUSort(uint ThreadIndex);
	void RunStats() const;
	uint GetQueryCount() const;

// Calibration
	void WriteCalibSample(FILE *f) const;
	void WriteCalibOutput(FILE *f) const;

public:
	virtual void OnSetup() {}
	virtual void OnAln(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA) {}
	virtual void OnAlnBA(uint ChainIndex1, uint ChainIndex2, DSSAligner &DA) {}

public:
	static void StaticThread(uint ThreadIndex, DBSearcher *ptrDBS);
	static void StaticThreadUSort(uint ThreadIndex, DBSearcher *ptrDBS);
	};
