#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"
#include "sort.h"

static vector<vector<byte> > s_Cols;
static vector<vector<byte> > s_Centroids;
static vector<uint> s_ClusterIdxs;
static uint s_K = UINT_MAX;

static void Load(const string &ChainsFN)
	{
	ProgressLog("Reading %s\n", ChainsFN.c_str());
	vector<PDBChain *> Chains;
	ReadChains(ChainsFN, Chains);
	const uint ChainCount = SIZE(Chains);
	DSS D;
	vector<vector<byte> > Profile;

	const uint FeatureCount = SIZE(DSSParams::m_Features);
	s_Cols.reserve(300*ChainCount);
	vector<byte> Col(FeatureCount);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Profiles");
		const PDBChain &Chain = *Chains[ChainIdx];
		D.Init(Chain);
		D.GetProfile(Profile);
		const uint L = Chain.GetSeqLength();
		asserta(L > 0);
		asserta(SIZE(Profile) == FeatureCount);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
				{
				const vector<byte> &Row = Profile[FeatureIdx];
				Col[FeatureIdx] = Row[Pos];
				}
			s_Cols.push_back(Col);
			}
		}
	ProgressLog("%u columns\n", SIZE(s_Cols));
	}

static void SetInitialCentroids()
	{
	set<uint> Idxs;
	const uint ColCount = SIZE(s_Cols);
	s_Centroids.clear();
	for (uint i = 0; i < s_K; ++i)
		{
		uint Idx = UINT_MAX;
		for (uint Try = 0; Try < 1000; ++Try)
			{
			Idx = randu32()%ColCount;
			if (Idxs.find(Idx) == Idxs.end())
				{
				Idxs.insert(Idx);
				s_Centroids.push_back(s_Cols[Idx]);
				break;
				}
			}
		asserta(Idx != UINT_MAX);
		}
	}

static float GetScoreColPair(const vector<byte> &Col1, const vector<byte> &Col2)
	{
	float Score = 0;
	const uint FeatureCount = SIZE(DSSParams::m_Features);
	assert(SIZE(Col1) == FeatureCount);
	assert(SIZE(Col2) == FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = DSSParams::m_Features[FeatureIdx];
		float **ScoreMx = DSSParams::m_ScoreMxs[F];
		byte ia = Col1[FeatureIdx];
		byte ib = Col2[FeatureIdx];
#if DEBUG
		uint AlphaSize = g_AlphaSizes2[F];
		assert(ia < AlphaSize && ib < AlphaSize);
#endif
		Score += ScoreMx[ia][ib];
		}
	return Score;
	}

static uint AssignCluster(uint ColIdx)
	{
	const vector<byte> &Col = s_Cols[ColIdx];
	asserta(SIZE(s_Centroids) == s_K);
	uint BestClusterIdx = UINT_MAX;
	float BestScore = -999;
	for (uint ClusterIdx = 0; ClusterIdx < s_K; ++ClusterIdx)
		{
		const vector<byte> &Center = s_Centroids[ClusterIdx];
		float Score = GetScoreColPair(Col, Center);
		if (Score > BestScore)
			{
			BestClusterIdx = ClusterIdx;
			BestScore = Score;
			}
		}
	asserta(BestClusterIdx != UINT_MAX);
	return BestClusterIdx;
	}

static uint AssignClusters()
	{
	const uint ColCount = SIZE(s_Cols);
	uint ChangeCount = 0;
	if (s_ClusterIdxs.empty())
		{
		s_ClusterIdxs.clear();
		s_ClusterIdxs.resize(ColCount, UINT_MAX);
		}
	else
		asserta(SIZE(s_ClusterIdxs) == ColCount);
	for (uint ColIdx = 0; ColIdx < ColCount; ++ColIdx)
		{
		uint ClusterIdx = AssignCluster(ColIdx);
		if (ClusterIdx != s_ClusterIdxs[ColIdx])
			{
			++ChangeCount;
			s_ClusterIdxs[ColIdx] = ClusterIdx;
			}
		}
	return ChangeCount;
	}

static void GetColIdxs(uint ClusterIdx, vector<uint> &ColIdxs)
	{
	ColIdxs.clear();
	const uint ColCount = SIZE(s_Cols);
	for (uint ColIdx = 0; ColIdx < ColCount; ++ColIdx)
		{
		if (s_ClusterIdxs[ColIdx] == ClusterIdx)
			ColIdxs.push_back(ColIdx);
		}
	}

static uint GetClusterSize(uint ClusterIdx)
	{
	const uint ColCount = SIZE(s_Cols);
	uint Size = 0;
	for (uint ColIdx = 0; ColIdx < ColCount; ++ColIdx)
		{
		if (s_ClusterIdxs[ColIdx] == ClusterIdx)
			++Size;
		}
	return Size;
	}

static byte GetModalLetter(const vector<uint> ColIdxs, uint FeatureIdx)
	{
	FEATURE F = DSSParams::m_Features[FeatureIdx];
	const uint AlphaSize = DSSParams::GetAlphaSize(F);
	vector<uint> Counts(AlphaSize, 0);
	const uint n = SIZE(ColIdxs);
	if (n == 0)
		return UINT8_MAX;
	for (uint i = 0; i < n; ++i)
		{
		uint ColIdx = ColIdxs[i];
		uint Letter = s_Cols[ColIdx][FeatureIdx];
		if (Letter < AlphaSize)
			Counts[Letter] += 1;
		else
			asserta(Letter == UINT_MAX);
		}
	uint BestLetter = UINT_MAX;
	uint BestCount = 0;
	for (uint Letter = 0; Letter < AlphaSize; ++Letter)
		{
		if (Counts[Letter] > BestCount)
			{
			BestCount = Counts[Letter];
			BestLetter = Letter;
			}
		}
	asserta(BestLetter < AlphaSize);
	return BestLetter;
	}

static void GetNewCentroid(uint ClusterIdx, vector<byte> &Centroid)
	{
	Centroid.clear();
	const uint FeatureCount = SIZE(DSSParams::m_Features);
	Centroid.resize(FeatureCount, UINT8_MAX);
	vector<uint> ColIdxs;
	GetColIdxs(ClusterIdx, ColIdxs);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		byte Letter = GetModalLetter(ColIdxs, FeatureIdx);
		Centroid[FeatureIdx] = Letter;
		}
	}

static void SetNewCentroids()
	{
	for (uint ClusterIdx = 0; ClusterIdx < s_K; ++ClusterIdx)
		GetNewCentroid(ClusterIdx, s_Centroids[ClusterIdx]);
	}

static void LogClusters()
	{
	const uint FeatureCount = SIZE(DSSParams::m_Features);
	Log("Cluster     Size");
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		FEATURE F = DSSParams::m_Features[FeatureIdx];
		Log("  %10.10s", FeatureToStr(F));
		}
	Log("\n");

	vector<uint> Sizes;
	for (uint ClusterIdx = 0; ClusterIdx < s_K; ++ClusterIdx)
		{
		uint Size = GetClusterSize(ClusterIdx);
		Sizes.push_back(Size);
		}
	vector<uint> Order(s_K);
	QuickSortOrderDesc(Sizes.data(), s_K, Order.data());
	uint SumSize = 0;
	for (uint j = 0; j < s_K; ++j)
		{
		uint ClusterIdx = Order[j];
		uint Size = Sizes[ClusterIdx];
		SumSize += Size;
		Log("%7u", ClusterIdx);
		Log("  %7u", Size);
		for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
			{
			byte Letter = s_Centroids[ClusterIdx][FeatureIdx];
			Log("  %10u", Letter);
			}
		Log("\n");
		}
	Log("%7.7s", "Total");
	Log("  %7u\n", SumSize);
	asserta(SumSize == SIZE(s_Cols));
	}

void cmd_nucluster()
	{
	const string &ChainsFN = g_Arg1;
	DSSParams::Init(DM_UseCommandLineOption);

	asserta(optset_k);
	s_K = opt(k);

	Load(ChainsFN);
	SetInitialCentroids();
	const uint ITERS = 100;
	for (uint Iter = 0; Iter < ITERS; ++Iter)
		{
		uint ChangeCount = AssignClusters();
		Progress("Iter %u, changes %u\n", Iter+1, ChangeCount);
		if (ChangeCount == 0)
			break;
		SetNewCentroids();
		}
	LogClusters();
	}
