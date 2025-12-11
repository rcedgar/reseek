#include "myutils.h"
#include "fragaligner.h"
#include "sort.h"
#include <set>

static uint s_M;
static uint s_AlphaSize;
static uint s_IdxCount;
static vector<PDBChain *> s_Chains;
static vector<uint *> s_DistsPtrVec;
static float *s_DistPairScoresPtr;

static vector<uint> s_ChainIdxs;
static vector<uint> s_FragStarts;

static vector<uint *> s_CentroidsDistsPtrVec;
static uint *s_ClusterIdxs;
static uint *s_PrevClusterIdxs;
static vector<vector<uint> > s_ClusterIdxToFragIdxs;

static mutex s_Lock;
static uint s_NextFragIdx;
static uint s_Chunk = 256;
static uint s_BandWidth = 2;
static uint s_DistN = 40;

static void SetRandomInitialCentroids()
	{
	s_CentroidsDistsPtrVec.clear();
	set<uint> DoneSet;
	const uint FragCount = SIZE(s_DistsPtrVec);
	for (uint i = 0; i < s_AlphaSize; ++i)
		{
		uint FragIdx = UINT_MAX;
		for (int Try = 0; Try < 100; ++Try)
			{
			FragIdx = randu32()%FragCount;
			if (DoneSet.find(FragIdx) != DoneSet.end())
				break;
			}
		asserta(FragIdx != UINT_MAX);
		DoneSet.insert(FragIdx);
		uint *DistsPtr = s_DistsPtrVec[FragIdx];
		s_CentroidsDistsPtrVec.push_back(DistsPtr);
		}
	}

static void LogCentroidDists()
	{
	Log("\n");
	Log("LogCentroidDists()\n");
	for (uint ClusterIdx = 0; ClusterIdx < s_AlphaSize; ++ClusterIdx)
		{
		Log("[%3u]  ", ClusterIdx);
		const uint *DistsPtr_Centroid = s_CentroidsDistsPtrVec[ClusterIdx];
		for (uint j = 0; j < s_IdxCount; ++j)
			Log(" %2u", DistsPtr_Centroid[j]);
		Log("\n");
		}
	}

static uint *GetRandomDistsPtr()
	{
	asserta(false);
	return 0;
	}

static void AssignMeanCentroid(uint ClusterIdx)
	{
	const vector<uint> &FragIdxs = s_ClusterIdxToFragIdxs[ClusterIdx];
	const uint n = SIZE(FragIdxs);
	if (n == 0)
		{
		s_CentroidsDistsPtrVec[ClusterIdx] = GetRandomDistsPtr();
		return;
		}

	vector<float> Sums(s_IdxCount);
	for (uint i = 0; i < n; ++i)
		{
		uint FragIdx = FragIdxs[i];
		const uint *DistsPtr = s_DistsPtrVec[FragIdx];
		for (uint j = 0; j < s_IdxCount; ++j)
			Sums[j] += DistsPtr[j];
		}

	for (uint j = 0; j < s_IdxCount; ++j)
		{
		uint Mean = uint(round(Sums[j]/n));
		s_CentroidsDistsPtrVec[ClusterIdx][j] = Mean;
		}
	}

static void AssignMeanCentroids()
	{
	for (uint ClusterIdx = 0; ClusterIdx < s_AlphaSize; ++ClusterIdx)
		AssignMeanCentroid(ClusterIdx);
	}

static uint AssignCluster(FragAligner &FA, uint FragIdx)
	{
	const uint *DistsPtr_Frag = s_DistsPtrVec[FragIdx];
	float BestScore = -999;
	uint BestClusterIdx = UINT_MAX;
	for (uint ClusterIdx = 0; ClusterIdx < s_AlphaSize; ++ClusterIdx)
		{
		const uint *DistsPtr_Centroid = s_CentroidsDistsPtrVec[ClusterIdx];
		float Score = FA.Align(DistsPtr_Frag, DistsPtr_Centroid);
		if (FragIdx == 0)
			{//@@
			Log("Frag 0 cluster %u score %.2f\n", ClusterIdx, Score);
			}//@@
		if (Score > BestScore)
			{
			BestScore = Score;
			BestClusterIdx = ClusterIdx;
			}
		}
	return BestClusterIdx;
	}

static void ThreadBody(uint ThreadIndex)
	{
	const uint FragCount = SIZE(s_DistsPtrVec);
	asserta(SIZE(s_CentroidsDistsPtrVec) == s_AlphaSize);
	FragAligner FA;
	FA.Init(s_M, s_BandWidth, s_DistN, s_DistPairScoresPtr);
	for (;;)
		{
		s_Lock.lock();
		uint NextFragIdx = s_NextFragIdx;
		if (s_NextFragIdx < FragCount)
			s_NextFragIdx += s_Chunk;
		s_Lock.unlock();
		if (s_NextFragIdx >= FragCount)
			break;

		uint IdxHi = min(NextFragIdx + s_Chunk, FragCount);
		vector<uint> ClusterIdxs;
		for (uint i = 0; i < s_Chunk; ++i)
			{
			uint FragIdx = i + NextFragIdx;
			if (FragIdx >= FragCount)
				break;
			uint ClusterIdx = AssignCluster(FA, FragIdx);
			ClusterIdxs.push_back(ClusterIdx);
			}

		s_Lock.lock();
		for (uint i = 0; i < s_Chunk; ++i)
			{
			uint FragIdx = i + NextFragIdx;
			if (FragIdx >= FragCount)
				break;
			uint ClusterIdx = ClusterIdxs[i];
			s_ClusterIdxs[FragIdx] = ClusterIdx;
			s_ClusterIdxToFragIdxs[ClusterIdx].push_back(FragIdx);
			}
		s_Lock.unlock();
		}
	}

static void AddFrag(FragAligner &FA, uint ChainIdx, uint Pos)
	{
	const PDBChain &Chain = *s_Chains[ChainIdx];
	uint *DistsPtr = FA.GetDistsPtr(Chain, Pos);
	s_ChainIdxs.push_back(ChainIdx);
	s_FragStarts.push_back(Pos);
	s_DistsPtrVec.push_back(DistsPtr);
	}

static void AddFrags(FragAligner &FA, uint ChainIdx)
	{
	const PDBChain &Chain = *s_Chains[ChainIdx];
	const uint L = Chain.GetSeqLength();
	for (uint Pos = randu32()%s_M; Pos + s_M <= L; Pos += s_M/4)
		AddFrag(FA, ChainIdx, Pos);
	}

static void LogClusterAssignsHead(uint n = 10)
	{
	const uint FragCount = SIZE(s_DistsPtrVec);
	Log("\n");
	Log("LogClusterAssignsHead()\n");
	for (uint FragIdx = 0; FragIdx < min(FragCount, n); ++FragIdx)
		Log("Frag %6u  cluster %3u\n", FragIdx, s_ClusterIdxs[FragIdx]);
	}

static void AssignClusters()
	{
	const uint FragCount = SIZE(s_DistsPtrVec);
	if (s_ClusterIdxs == 0)
		{
		s_ClusterIdxs = myalloc(uint, FragCount);
		s_PrevClusterIdxs = myalloc(uint, FragCount);
		}
	s_ClusterIdxToFragIdxs.clear();
	s_ClusterIdxToFragIdxs.resize(s_AlphaSize);

	vector<thread *> ts;
	uint ThreadCount = GetRequestedThreadCount();
	s_NextFragIdx = 0;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	}

void cmd_build_s3e_library()
	{
	const string &ChainsFN = g_Arg1;
	asserta(optset_subsample);
	s_M = opt(m);

	asserta(optset_alpha_size);
	s_AlphaSize = opt(alpha_size);

	asserta(optset_m);
	s_M = opt(m);

	s_BandWidth = 2;
	s_DistN = 20;
	s_DistPairScoresPtr = FragAligner::MakeDistPairScoresPtr(s_DistN);

	FragAligner FA;
	FA.Init(s_M, s_BandWidth, s_DistN, s_DistPairScoresPtr);
	s_IdxCount = FA.m_IdxCount;

	Progress("Reading chains...");
	ReadChains(ChainsFN, s_Chains);
	const uint ChainCount = SIZE(s_Chains);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Add frags");
		AddFrags(FA, ChainIdx);
		}
	uint FragCount = SIZE(s_ChainIdxs);
	asserta(SIZE(s_FragStarts) == FragCount);
	asserta(SIZE(s_DistsPtrVec) == FragCount);
	float FragsPerChain = float(FragCount)/ChainCount;
	ProgressLog("%u frags, %.2f frags/chain\n", FragCount, FragsPerChain);

	SetRandomInitialCentroids();
	AssignClusters();
	uint BestChangeCount = FragCount;
	uint BestIter = UINT_MAX;
	for (uint Iter = 0; Iter < 1000; ++Iter)
		{
		//LogCentroidDists();
		//LogClusterAssignsHead();
		memcpy(s_PrevClusterIdxs, s_ClusterIdxs, FragCount*sizeof(uint));
		AssignMeanCentroids();
		AssignClusters();
		uint ChangeCount = 0;
		for (uint FragIdx = 0; FragIdx < FragCount; ++FragIdx)
			{
			if (s_PrevClusterIdxs[FragIdx] != s_ClusterIdxs[FragIdx])
				++ChangeCount;
			}
		if (ChangeCount < BestChangeCount)
			{
			BestChangeCount = ChangeCount;
			BestIter = Iter;
			}
		ProgressLog("[%u] changes %u (%u, %u)\n",
			Iter+1, ChangeCount, BestChangeCount, BestIter+1);
		if (ChangeCount == 0 || Iter - BestIter > 1000)
			{
			ProgressLog("Converged\n");
			break;
			}
		}
	ProgressLog("Done\n");
	}
