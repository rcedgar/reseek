#include "myutils.h"
#include "pdbchain.h"
#include "dss.h"
#include "sort.h"
#include "binner.h"
#include "trainer.h"
#include <omp.h>

static vector<vector<byte> > s_Cols;
static vector<vector<byte> > s_Centroids;
static vector<uint> s_ClusterIdxs;
static vector<uint> s_ClusterSizes;
static uint s_K = UINT_MAX;

static vector<vector<byte> > s_Profile_Q;
static vector<vector<byte> > s_Profile_T;

static void GetProfileCol(const vector<vector<byte> > &Profile,
	uint Pos, vector<byte> &Col)
	{
	Col.clear();
	const uint FeatureCount = SIZE(DSSParams::m_Features);
	Col.reserve(FeatureCount);
	asserta(SIZE(Profile) == FeatureCount);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		const vector<byte> &Row = Profile[FeatureIdx];
		Col.push_back(Row[Pos]);
		}
	}

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
	vector<byte> Col;
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
			GetProfileCol(Profile, Pos, Col);
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

static uint AssignCluster(const vector<byte> &Col)
	{
	asserta(SIZE(s_Centroids) == s_K);
	asserta(SIZE(s_ClusterSizes) == s_K);
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

		uint ColCount = SIZE(s_Cols);
		uint IdealSize = ColCount/s_K;
		s_ClusterSizes.clear();
		s_ClusterSizes.resize(s_K, IdealSize);
		}
	else
		asserta(SIZE(s_ClusterIdxs) == ColCount);
	for (uint ColIdx = 0; ColIdx < ColCount; ++ColIdx)
		{
		const vector<byte> &Col = s_Cols[ColIdx];
		uint ClusterIdx = AssignCluster(Col);
		if (ClusterIdx != s_ClusterIdxs[ColIdx])
			{
			++ChangeCount;
			s_ClusterIdxs[ColIdx] = ClusterIdx;
			s_ClusterSizes[ClusterIdx] += 1;
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

static void GetCentroid_ModalLetter(
	const vector<uint> &ColIdxs, vector<byte> &Centroid)
	{
	Centroid.clear();
	const uint FeatureCount = SIZE(DSSParams::m_Features);
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		{
		byte Letter = GetModalLetter(ColIdxs, FeatureIdx);
		Centroid.push_back(Letter);
		}
	}

static float GetAvgScore(const vector<byte> &Col,
	const vector<uint> &ColIdxs)
	{
	uint n = SIZE(ColIdxs);
	if (n == 0)
		return 0;
	float Total = 0;
	for (uint i = 0; i < n; ++i)
		Total += GetScoreColPair(Col, s_Cols[ColIdxs[i]]);
	return Total/n;
	}

static void GetCentroid_MaximizeScore(
	const vector<uint> &ColIdxs, vector<byte> &Centroid)
	{
	const uint FeatureCount = SIZE(DSSParams::m_Features);
	Centroid.clear();
	Centroid.resize(FeatureCount, 0);
	vector<vector<byte> > Cols;
	const uint N = SIZE(ColIdxs);
	if (N == 0)
		return;
	float BestAvgScore = -999;
#pragma omp parallel for num_threads(32)
	for (int i = 0; i < int(N); ++i)
		{
		const vector<byte> &Col = s_Cols[ColIdxs[randu32()%N]];
		float Score = GetAvgScore(Col, ColIdxs);
#pragma omp critical
		{
		if (Score > BestAvgScore)
			{
			BestAvgScore = Score;
			Centroid = Col;
			}
		}
		}
	}

static void GetNewCentroid(uint ClusterIdx, vector<byte> &Centroid)
	{
	Centroid.clear();
	const uint FeatureCount = SIZE(DSSParams::m_Features);
	Centroid.resize(FeatureCount, UINT8_MAX);
	vector<uint> ColIdxs;
	GetColIdxs(ClusterIdx, ColIdxs);
	if (opt(modal))
		return GetCentroid_ModalLetter(ColIdxs, Centroid);
	return GetCentroid_MaximizeScore(ColIdxs, Centroid);
	}

static void SetNewCentroids()
	{
	for (uint ClusterIdx = 0; ClusterIdx < s_K; ++ClusterIdx)
		{
		ProgressStep(ClusterIdx, s_K, "New centroids");
		GetNewCentroid(ClusterIdx, s_Centroids[ClusterIdx]);
		}
	}

static void LogClusters()
	{
	const uint FeatureCount = SIZE(DSSParams::m_Features);
	Log("\n");
	Log("Cluster     Size    AvgScore   SelfScore");
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
	float TotalAvgScore = 0;
	for (uint j = 0; j < s_K; ++j)
		{
		uint ClusterIdx = Order[j];
		uint Size = Sizes[ClusterIdx];
		SumSize += Size;
		vector<uint> ColIdxs;
		GetColIdxs(ClusterIdx, ColIdxs);
		const vector<byte> &Centroid = s_Centroids[ClusterIdx];
		float AvgScore = GetAvgScore(Centroid, ColIdxs);
		TotalAvgScore += AvgScore;
		float SelfScore = GetScoreColPair(Centroid, Centroid);
		Log("%7u", ClusterIdx);
		Log("  %7u", Size);
		Log("  %10.3g", AvgScore);
		Log("  %10.3g", SelfScore);
		for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
			{
			byte Letter = s_Centroids[ClusterIdx][FeatureIdx];
			Log("  %10u", Letter);
			}
		Log("\n");
		}
	Log("%7.7s", "Total");
	Log("  %7u", SumSize);
	Log("  %10.3g\n", TotalAvgScore/s_K);
	asserta(SumSize == SIZE(s_Cols));
	}

static void ClustersToTsv(FILE *f)
	{
	if (f == 0)
		return;
	vector<uint> Sizes;
	for (uint ClusterIdx = 0; ClusterIdx < s_K; ++ClusterIdx)
		{
		uint Size = GetClusterSize(ClusterIdx);
		Sizes.push_back(Size);
		}
	vector<uint> Order(s_K);
	QuickSortOrderDesc(Sizes.data(), s_K, Order.data());
	uint SumSize = 0;
	float TotalAvgScore = 0;

	const uint FeatureCount = SIZE(DSSParams::m_Features);
	fprintf(f, "AS\t%u\n", s_K);
	fprintf(f, "Cluster");
	fprintf(f, "\tSize");
	fprintf(f, "\tAvgScore");
	fprintf(f, "\tSelfScore");
	for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
		fprintf(f, "\t%s", FeatureToStr(DSSParams::m_Features[FeatureIdx]));
	fprintf(f, "\n");
	for (uint j = 0; j < s_K; ++j)
		{
		uint ClusterIdx = Order[j];
		uint Size = Sizes[ClusterIdx];
		SumSize += Size;
		vector<uint> ColIdxs;
		GetColIdxs(ClusterIdx, ColIdxs);
		const vector<byte> &Centroid = s_Centroids[ClusterIdx];
		float AvgScore = GetAvgScore(Centroid, ColIdxs);
		TotalAvgScore += AvgScore;
		float SelfScore = GetScoreColPair(Centroid, Centroid);

		fprintf(f, "%u", j);
		fprintf(f, "\t%u", Size);
		fprintf(f, "\t%.3g", AvgScore);
		fprintf(f, "\t%.3g", SelfScore);
		for (uint FeatureIdx = 0; FeatureIdx < FeatureCount; ++FeatureIdx)
			{
			byte Letter = s_Centroids[ClusterIdx][FeatureIdx];
			fprintf(f, "\t%u", Letter);
			}
		fprintf(f, "\n");
		}
	asserta(SumSize == SIZE(s_Cols));
	}

static void LogScoreBins()
	{
	vector<float> Scores;
	const uint N = 1000000;
	const uint ColCount = SIZE(s_Cols);
	Scores.reserve(N);
	for (uint i = 0; i < N; ++i)
		{
		uint Idx1 = randu32()%ColCount;
		uint Idx2 = randu32()%ColCount;
		const vector<byte> &Col1 = s_Cols[Idx1];
		const vector<byte> &Col2 = s_Cols[Idx2];
		float Score = GetScoreColPair(Col1, Col2);
		Scores.push_back(Score);
		}
	Binner<float> B(Scores, 32);
	Log("\n");
	B.ToHist(g_fLog);
	}

static void TrainerOnPair(
  const Trainer &T, uint ChainIdxQ, uint ChainIdxT,
  const vector<uint> &PosQs, const vector<uint> &PosRs)
	{
	DSS D;
	const PDBChain &ChainQ = T.GetChain(ChainIdxQ);
	const PDBChain &ChainT = T.GetChain(ChainIdxT);
	D.Init(ChainQ);
	D.GetProfile(s_Profile_Q);

	D.Init(ChainT);
	D.GetProfile(s_Profile_T);
	}

static void ChainToFasta(FILE *f, const PDBChain &Chain)
	{
	if (f == 0)
		return;
	DSS D;
	D.Init(Chain);
	vector<vector<byte> > Profile;
	D.GetProfile(Profile);
	const uint L = Chain.GetSeqLength();
	vector<byte> Col;
	string HexSeq;
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		GetProfileCol(Profile, Pos, Col);
		uint ClusterIdx = AssignCluster(Col);
		asserta(ClusterIdx < s_K);
		Psa(HexSeq, "%02x", ClusterIdx);
		}
	SeqToFasta(f, Chain.m_Label, HexSeq);
	}

static void ChainsToFasta(const vector<PDBChain *> &Chains, const string &FN)
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);
	const uint ChainCount = SIZE(Chains);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "%s", FN.c_str());
		const PDBChain &Chain = *Chains[ChainIdx];
		ChainToFasta(f, Chain);
		}
	CloseStdioFile(f);
	}

static void TrainerAlphaCol(
  const Trainer &T, uint PosQ, uint PosT,
  uint &LetterQ, uint &LetterT)
	{
	const uint FeatureCount = SIZE(DSSParams::m_Features);
	asserta(SIZE(s_Profile_Q) == FeatureCount);
	asserta(SIZE(s_Profile_T) == FeatureCount);

	vector<byte> Col_T, Col_Q;
	GetProfileCol(s_Profile_T, PosT, Col_T);
	GetProfileCol(s_Profile_Q, PosQ, Col_Q);

	LetterQ = AssignCluster(Col_Q);
	LetterT = AssignCluster(Col_T);
	}

void cmd_nucluster()
	{
	const string &ChainsFN = g_Arg1;
	DSSParams::Init(DM_UseCommandLineOption);

	asserta(optset_k);
	s_K = opt(k);

	Load(ChainsFN);
	LogScoreBins();
	SetInitialCentroids();
	const uint ITERS = 100;
	const uint IMPROVED_WINDOW = 10;
	vector<uint> ChangeCounts;
	uint BestIter = UINT_MAX;
	uint MinChangeCount = UINT_MAX;
	vector<vector<byte> > BestCentroids;
	for (uint Iter = 0; Iter < ITERS; ++Iter)
		{
		uint ChangeCount = AssignClusters();
		if (ChangeCount == 0)
			{
			Progress("Iter %u, converged no changes\n", Iter+1);
			break;
			}
		ChangeCounts.push_back(ChangeCount);
		int dc = int(MinChangeCount) - int(ChangeCount);
		int di = int(Iter) - int(BestIter);
		Progress("Iter %u, changes %u [%+d, %u]\n", Iter+1, ChangeCount, dc, di);
		if (ChangeCount < MinChangeCount)
			{
			BestIter = Iter;
			MinChangeCount = ChangeCount;
			BestCentroids = s_Centroids;
			}
		if (Iter - BestIter > IMPROVED_WINDOW)
			{
			ProgressLog("Converged, best iter %u\n", BestIter+1);
			s_Centroids = BestCentroids;
			break;
			}
		SetNewCentroids();
		}
	LogClusters();

	if (optset_db && optset_traintps)
		{
		asserta(optset_output);
		Progress("Training log-odds score mx\n");
		Trainer Tr;
		Tr.Init(opt(traintps), opt(db));
		LogOdds LO;
		LO.m_UseUnalignedBackground = false;
		Tr.TrainLogOdds(s_K, TrainerOnPair, TrainerAlphaCol, LO);
		uint pc = 3;
		if (optset_psuedocount)
			pc = opt(psuedocount);
		LO.AddPseudoCount(pc);

		vector<vector<float> > ScoreMx;
		double ExpectedScore = LO.GetLogOddsMx(ScoreMx);
		LO.ToTsv(g_fLog);
		LO.MxToSrc(g_fLog, "Nu", ScoreMx);
		ProgressLog("ES %.3g\n", ExpectedScore);

		FILE *fOut = CreateStdioFile(opt(output));
		string ParamStr;
		DSSParams::GetParamsStr(ParamStr);
		string CmdLine;
		GetCmdLine(CmdLine);

		fprintf(fOut, "# %s\n", CmdLine.c_str());
		fprintf(fOut, "params\t%s\n", opt(params));
		fprintf(fOut, "params2\t%s\n", ParamStr.c_str());
		fprintf(fOut, "pseudocount\t%u\n", pc);
		LO.MxToTsv(fOut, "Nu", ScoreMx);
		ClustersToTsv(fOut);
		CloseStdioFile(fOut);
		if (optset_fasta)
			ChainsToFasta(Tr.m_Chains, opt(fasta));
		}
	}
