#pragma once

#include "fragaligner.h"

class PDBChain;

class SSSLib
	{
public:
	uint m_M = UINT_MAX;
	uint m_AlphaSize = UINT_MAX;
	uint m_IdxCount = UINT_MAX;
	uint m_NextFragIdx = 0;
	uint m_Chunk = 256;
	uint m_BandWidth = 2;
	uint m_DistN = 40;

	vector<PDBChain *> m_Chains;

// Per training fragment vectors
	vector<uint *> m_DistsPtrVec;
	vector<uint> m_FragChainIdxs;
	vector<uint> m_FragStarts;
	uint *m_ClusterIdxs = 0;
	uint *m_PrevClusterIdxs = 0;

// Per-centroid vectors
	vector<uint *> m_CentroidsDistsPtrVec;
	vector<vector<uint> > m_ClusterIdxToFragIdxs;
	vector<uint> m_RepFragIdxs;

// Scoring matrix for pair of integer distances
	float *m_DistPairScoresPtr = 0;

	mutex m_Lock;
	FragAligner m_FA;

	vector<vector<byte> > m_LettersVec;

public:
	void SetRandomInitialCentroids();
	void LogCentroidDists();
	void AssignMeanCentroid(uint ClusterIdx);
	uint *GetRandomDistsPtr();
	void AssignMeanCentroids();
	void AssignClusters();
	uint AssignCluster(FragAligner &FA, uint FragIdx);
	void AddFrag(FragAligner &FA, uint ChainIdx, uint Pos);
	void AddFrags(FragAligner &FA, uint ChainIdx);
	void LogClusterAssignsHead(uint n = 10);
	void ThreadBody(uint ThreadIndex);
	void SetParams(uint AlphaSize, uint M, uint BandWidth, uint DistN);
	void SetTrainingFrags();
	void Train();
	void LoadChains(const string &FN);
	void AssertSame(FragAligner &FA, uint FragIdx, uint ClusterIdx);
	void AssertSames();
	uint GetFragCount() const { return SIZE(m_DistsPtrVec); }
	void AssignReps();
	void GetClusterSizes(vector<uint> &Sizes, vector<uint> &Order) const;
	void LogClusters() const;
	void AssignLetters();
	byte GetLetter(const uint *DistsPtr, uint L, uint Start);
	void ToFasta(const string &FN) const;
	char LetterToChar(byte Letter) const;

public:
	static void StaticThreadBody(SSSLib *Lib, uint ThreadIndex);
	};
