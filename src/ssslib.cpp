#include "myutils.h"
#include "ssslib.h"
#include "pdbchain.h"
#include "fragaligner.h"
#include "sort.h"
#include "alpha.h"

void SSSLib::SetRandomInitialCentroids()
	{
	m_CentroidsDistsPtrVec.clear();
	set<uint> DoneSet;
	const uint FragCount = SIZE(m_DistsPtrVec);
	for (uint i = 0; i < m_AlphaSize; ++i)
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
		uint *DistsPtr = m_DistsPtrVec[FragIdx];
		m_CentroidsDistsPtrVec.push_back(DistsPtr);
		}
	}

void SSSLib::LogCentroidDists()
	{
	Log("\n");
	Log("LogCentroidDists()\n");
	for (uint ClusterIdx = 0; ClusterIdx < m_AlphaSize; ++ClusterIdx)
		{
		Log("[%3u]  ", ClusterIdx);
		const uint *DistsPtr_Centroid = m_CentroidsDistsPtrVec[ClusterIdx];
		for (uint j = 0; j < m_IdxCount; ++j)
			Log(" %2u", DistsPtr_Centroid[j]);
		Log("\n");
		}
	}

uint *SSSLib::GetRandomDistsPtr()
	{
	uint *DistsPtr = myalloc(uint, m_IdxCount);
	for (uint i = 0; i < m_IdxCount; ++i)
		{
		uint ClusterIdx = randu32()%m_AlphaSize;
		DistsPtr[i] = m_CentroidsDistsPtrVec[ClusterIdx][i];
		}
	return DistsPtr;
	}

void SSSLib::AssignMeanCentroid(uint ClusterIdx)
	{
	const vector<uint> &FragIdxs = m_ClusterIdxToFragIdxs[ClusterIdx];
	const uint n = SIZE(FragIdxs);
	if (n == 0)
		{
		m_CentroidsDistsPtrVec[ClusterIdx] = GetRandomDistsPtr();
		return;
		}

	vector<float> Sums(m_IdxCount);
	for (uint i = 0; i < n; ++i)
		{
		uint FragIdx = FragIdxs[i];
		const uint *DistsPtr = m_DistsPtrVec[FragIdx];
		for (uint j = 0; j < m_IdxCount; ++j)
			Sums[j] += DistsPtr[j];
		}

	for (uint j = 0; j < m_IdxCount; ++j)
		{
		uint Mean = uint(round(Sums[j]/n));
		m_CentroidsDistsPtrVec[ClusterIdx][j] = Mean;
		}
	}

void SSSLib::AssignMeanCentroids()
	{
	for (uint ClusterIdx = 0; ClusterIdx < m_AlphaSize; ++ClusterIdx)
		AssignMeanCentroid(ClusterIdx);
	}

uint SSSLib::AssignCluster(FragAligner &FA, uint FragIdx)
	{
	const uint *DistsPtr_Frag = m_DistsPtrVec[FragIdx];
	float BestScore = -999;
	uint BestClusterIdx = UINT_MAX;
	for (uint ClusterIdx = 0; ClusterIdx < m_AlphaSize; ++ClusterIdx)
		{
		const uint *DistsPtr_Centroid = m_CentroidsDistsPtrVec[ClusterIdx];
		float Score = FA.Align(DistsPtr_Frag, DistsPtr_Centroid);
		if (Score > BestScore)
			{
			BestScore = Score;
			BestClusterIdx = ClusterIdx;
			}
		}
	return BestClusterIdx;
	}

void SSSLib::AssertSame(FragAligner &FA, uint FragIdx1, uint FragIdx2)
	{
	uint ChainIdx1 = m_FragChainIdxs[FragIdx1];
	uint ChainIdx2 = m_FragChainIdxs[FragIdx2];

	uint Start1 = m_FragStarts[FragIdx1];
	uint Start2 = m_FragStarts[FragIdx2];

	const PDBChain &Chain1 = *m_Chains[ChainIdx1];
	const PDBChain &Chain2 = *m_Chains[ChainIdx2];

	const uint L1 = Chain1.GetSeqLength();
	const uint L2 = Chain2.GetSeqLength();

	const uint *DistsPtr_Frag1 = m_DistsPtrVec[FragIdx1];
	const uint *DistsPtr_Frag2 = m_DistsPtrVec[FragIdx2];

	const uint *DistsPtr_Chain1 = FA.GetDistsPtrChain(Chain1);
	const uint *DistsPtr_Chain2 = FA.GetDistsPtrChain(Chain2);

	float ScoreSub = FA.AlignSubchain(DistsPtr_Chain1, L1, Start1, DistsPtr_Frag2);
	float ScoreFrag = FA.Align(DistsPtr_Frag1, DistsPtr_Frag2);
	asserta(feq(ScoreSub, ScoreFrag));
	}

void SSSLib::AssertSames()
	{
	uint FragCount = GetFragCount();
	for (uint Iter = 0; Iter < 1000; ++Iter)
		{
		uint i = randu32()%FragCount;
		uint j = randu32()%FragCount;
		AssertSame(m_FA, i, j);
		}
	ProgressLog("AssertSame ok\n");
	}

void SSSLib::ThreadBody(uint ThreadIndex)
	{
	const uint FragCount = SIZE(m_DistsPtrVec);
	asserta(SIZE(m_CentroidsDistsPtrVec) == m_AlphaSize);
	FragAligner FA;
	FA.Init(m_FragL, m_BandWidth, m_DistN, m_DistPairScoresPtr);
	for (;;)
		{
		m_Lock.lock();
		uint NextFragIdx = m_NextFragIdx;
		if (m_NextFragIdx < FragCount)
			m_NextFragIdx += m_ThreadBatch;
		m_Lock.unlock();
		if (m_NextFragIdx >= FragCount)
			break;

		uint IdxHi = min(NextFragIdx + m_ThreadBatch, FragCount);
		vector<uint> ClusterIdxs;
		for (uint i = 0; i < m_ThreadBatch; ++i)
			{
			uint FragIdx = i + NextFragIdx;
			if (FragIdx >= FragCount)
				break;
			uint ClusterIdx = AssignCluster(FA, FragIdx);
			ClusterIdxs.push_back(ClusterIdx);
			}

		m_Lock.lock();
		for (uint i = 0; i < m_ThreadBatch; ++i)
			{
			uint FragIdx = i + NextFragIdx;
			if (FragIdx >= FragCount)
				break;
			uint ClusterIdx = ClusterIdxs[i];
			m_ClusterIdxs[FragIdx] = ClusterIdx;
			m_ClusterIdxToFragIdxs[ClusterIdx].push_back(FragIdx);
			}
		m_Lock.unlock();
		}
	}

void SSSLib::StaticThreadBody(SSSLib *Lib, uint ThreadIndex)
	{
	Lib->ThreadBody(ThreadIndex);
	}

void SSSLib::AddFrag(FragAligner &FA, uint ChainIdx, uint Pos)
	{
	const PDBChain &Chain = *m_Chains[ChainIdx];
	uint *DistsPtr = FA.GetDistsPtrFrag(Chain, Pos);
	m_FragChainIdxs.push_back(ChainIdx);
	m_FragStarts.push_back(Pos);
	m_DistsPtrVec.push_back(DistsPtr);
	}

void SSSLib::AddFrags(FragAligner &FA, uint ChainIdx)
	{
	const PDBChain &Chain = *m_Chains[ChainIdx];
	const uint L = Chain.GetSeqLength();
	for (uint Pos = randu32()%m_FragL; Pos + m_FragL <= L; Pos += m_FragStep)
		AddFrag(FA, ChainIdx, Pos);
	}

void SSSLib::LogClusterAssignsHead(uint n)
	{
	const uint FragCount = SIZE(m_DistsPtrVec);
	Log("\n");
	Log("LogClusterAssignsHead()\n");
	for (uint FragIdx = 0; FragIdx < min(FragCount, n); ++FragIdx)
		Log("Frag %6u  cluster %3u\n", FragIdx, m_ClusterIdxs[FragIdx]);
	}

void SSSLib::AssignReps()
	{
	m_RepFragIdxs.clear();
	const uint FragCount = GetFragCount();
	for (uint ClusterIdx = 0; ClusterIdx < m_AlphaSize; ++ClusterIdx)
		{
		float BestScore = -999;
		uint BestFragIdx = UINT_MAX;
		const uint *DistsPtr_Centroid = m_CentroidsDistsPtrVec[ClusterIdx];
		for (uint FragIdx = 0; FragIdx < FragCount; ++FragIdx)
			{
			const uint *DistsPtr_Frag = m_DistsPtrVec[FragIdx];
			float Score = m_FA.Align(DistsPtr_Frag, DistsPtr_Centroid);
			if (Score > BestScore)
				{
				BestScore = Score;
				BestFragIdx = FragIdx;
				}
			}
		m_RepFragIdxs.push_back(BestFragIdx);
		}
	}

void SSSLib::AssignClusters()
	{
	const uint FragCount = SIZE(m_DistsPtrVec);
	if (m_ClusterIdxs == 0)
		{
		m_ClusterIdxs = myalloc(uint, FragCount);
		m_PrevClusterIdxs = myalloc(uint, FragCount);
		}
	m_ClusterIdxToFragIdxs.clear();
	m_ClusterIdxToFragIdxs.resize(m_AlphaSize);

	vector<thread *> ts;
	uint ThreadCount = GetRequestedThreadCount();
	m_NextFragIdx = 0;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadBody, this, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	}

void SSSLib::SetTrainingFrags()
	{
	const uint ChainCount = SIZE(m_Chains);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Add frags");
		AddFrags(m_FA, ChainIdx);
		}
	uint FragCount = SIZE(m_FragChainIdxs);
	asserta(SIZE(m_FragStarts) == FragCount);
	asserta(SIZE(m_DistsPtrVec) == FragCount);
	float FragsPerChain = float(FragCount)/ChainCount;
	ProgressLog("%u frags, %.2f frags/chain\n", FragCount, FragsPerChain);
	}

void SSSLib::SetParams(uint AlphaSize, uint M, uint BandWidth,
	uint DistN, uint FragStep)
	{
	m_AlphaSize = AlphaSize;
	m_FragL = M;
	m_BandWidth = BandWidth;
	m_DistN = DistN;
	m_FragStep = FragStep;

	m_DistPairScoresPtr = FragAligner::MakeDistPairScoresPtr(m_DistN);

	m_FA.Init(M, BandWidth, DistN, m_DistPairScoresPtr);
	m_IdxCount = m_FA.m_IdxCount;
	}

void SSSLib::Train()
	{
	const uint FragCount = SIZE(m_FragChainIdxs);
	SetRandomInitialCentroids();
	AssignClusters();
	uint BestChangeCount = FragCount;
	uint BestIter = UINT_MAX;
	for (uint Iter = 0; Iter < 1000; ++Iter)
		{
		memcpy(m_PrevClusterIdxs, m_ClusterIdxs, FragCount*sizeof(uint));
		AssignMeanCentroids();
		AssignClusters();
		uint ChangeCount = 0;
		for (uint FragIdx = 0; FragIdx < FragCount; ++FragIdx)
			{
			if (m_PrevClusterIdxs[FragIdx] != m_ClusterIdxs[FragIdx])
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

void SSSLib::GetClusterSizes(vector<uint> &Sizes, vector<uint> &Order) const
	{
	for (uint i = 0; i < m_AlphaSize; ++i)
		Sizes.push_back(SIZE(m_ClusterIdxToFragIdxs[i]));

	Order.resize(m_AlphaSize);
	QuickSortOrderDesc(Sizes.data(), m_AlphaSize, Order.data());
	}

void SSSLib::LoadChains(const string &FN)
	{
	ReadChains(FN, m_Chains);
	}

void SSSLib::LogClusters() const
	{
	vector<uint> Sizes, Order;
	GetClusterSizes(Sizes, Order);

	uint FragCount = GetFragCount();
	for (uint k = 0; k < m_AlphaSize; ++k)
		{
		uint ClusterIdx = Order[k];
		uint Size = Sizes[ClusterIdx];
		double Pct = GetPct(Size, FragCount);
		uint FragIdx = m_RepFragIdxs[ClusterIdx];
		uint ChainIdx = m_FragChainIdxs[FragIdx];
		uint Start = m_FragStarts[FragIdx];
		const string &Label = m_Chains[ChainIdx]->m_Label;
		Log("[%3u]", k);
		Log("  %7u", Size);
		Log("  %4.1f", Pct);
		Log("  %s(%u)", Label.c_str(), Start);
		Log("\n");
		}
	}

uint SSSLib::GetMinClusterSize() const
	{
	uint MinSize = UINT_MAX;
	for (uint i = 0; i < m_AlphaSize; ++i)
		MinSize = min(MinSize, SIZE(m_ClusterIdxToFragIdxs[i]));
	return MinSize;
	}

byte SSSLib::GetLetter(const uint *DistsPtr, uint L, uint Pos)
	{
	int iStart = int(Pos) - int(m_FragL)/2;
	if (iStart < 0)
		return 0xff;
	int iEnd = iStart + int(m_FragL) - 1;
	if (iEnd >= int(L))
		return 0xff;
	uint Start = uint(iStart);

	float BestScore = -999;
	byte BestLetter = 0xff;
	for (uint Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		const uint *DistsPtr_Centroid = m_CentroidsDistsPtrVec[Letter];
		float Score = m_FA.AlignSubchain(DistsPtr, L, Start, DistsPtr_Centroid);
		if (Score > BestScore)
			{
			BestLetter = Letter;
			BestScore = Score;
			}
		}
	return BestLetter;
	}

void SSSLib::AssignLetters()
	{
	const uint ChainCount = SIZE(m_Chains);
	m_LettersVec.clear();
	m_LettersVec.resize(ChainCount);
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Assign letters");
		const PDBChain &Chain = *m_Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		vector<byte> &Letters = m_LettersVec[ChainIdx];
		Letters.reserve(L);
		uint *DistsPtr = m_FA.GetDistsPtrChain(Chain);
		for (uint Pos = 0; Pos < L; ++Pos)
			Letters.push_back(GetLetter(DistsPtr, L, Pos));
		myfree(DistsPtr);
		}
	}

char SSSLib::LetterToChar(byte Letter) const
	{
	if (Letter == 0xff)
		return '*';
	asserta(Letter < 36);
	return g_LetterToCharMu[Letter];
	}

void SSSLib::ToSpec(const string &FN) const
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);
	vector<uint> Sizes, Order;
	GetClusterSizes(Sizes, Order);
	string CmdLine;
	GetCmdLine(CmdLine);
	fprintf(f, "# %s\n", CmdLine.c_str());
	fprintf(f, "alpha_size=%u\n", m_AlphaSize);
	fprintf(f, "fragl=%u\n", m_FragL);
	fprintf(f, "bandwidth=%u\n", m_BandWidth);
	fprintf(f, "fragstep=%u\n", m_FragStep);
	fprintf(f, "distn=%u\n", m_DistN);
	fprintf(f, "idxcount=%u\n", m_IdxCount);
	for (uint k = 0; k < m_AlphaSize; ++k)
		{
		uint ClusterIdx = Order[k];
		uint Size = Sizes[ClusterIdx];
		uint FragIdx = m_RepFragIdxs[ClusterIdx];
		uint ChainIdx = m_FragChainIdxs[FragIdx];
		uint Start = m_FragStarts[FragIdx];
		const string &Label = m_Chains[ChainIdx]->m_Label;
		fprintf(f, "%u", k);
		fprintf(f, "\t%7u", Size);
		fprintf(f, "\t%s(%u)", Label.c_str(), Start);
		const uint *DistsPtr = m_CentroidsDistsPtrVec[ClusterIdx];
		for (uint Idx = 0; Idx < m_IdxCount; ++Idx)
			fprintf(f, "\t%u", DistsPtr[Idx]);
		fprintf(f, "\n");
		}
	CloseStdioFile(f);
	}

void SSSLib::ToFasta(const string &FN) const
	{
	if (FN == "")
		return;
	FILE *f = CreateStdioFile(FN);
	const uint ChainCount = SIZE(m_Chains);
	asserta(SIZE(m_LettersVec) == ChainCount);
	string Tmp;
	for (uint ChainIdx = 0; ChainIdx < ChainCount; ++ChainIdx)
		{
		ProgressStep(ChainIdx, ChainCount, "Writing %s", FN.c_str());
		const PDBChain &Chain = *m_Chains[ChainIdx];
		const uint L = Chain.GetSeqLength();
		const string &Label = Chain.m_Label;
		const vector<byte> &Letters = m_LettersVec[ChainIdx];
		asserta(SIZE(Letters) == L);
		string Seq;
		if (m_AlphaSize <= 36)
			{
			for (uint Pos = 0; Pos < L; ++Pos)
				Seq += LetterToChar(Letters[Pos]);
			}
		else
			{
			for (uint Pos = 0; Pos < L; ++Pos)
				{
				Ps(Tmp, "%02x", Letters[Pos]);
				Seq += Tmp;
				}
			}
		SeqToFasta(f, Label, Seq);
		}
	CloseStdioFile(f);
	}

void cmd_cluster_sss()
	{
	const string &ChainsFN = g_Arg1;
	asserta(optset_alpha_size);
	asserta(optset_fragl);

	uint AlphaSize = opt(alpha_size);
	uint FragL = opt(fragl);
	uint BandWidth = 2;
	uint DistN = 20;
	uint FragStep = 8;
	uint MinSize = 100;

	if (optset_bandwidth)
		BandWidth = opt(bandwidth);
	if (optset_bandwidth)
		FragStep = opt(fragstep);
	if (optset_distn)
		DistN = opt(distn);
	if (optset_minsize)
		DistN = opt(minsize);

	SSSLib Lib;
	Lib.SetParams(AlphaSize, FragL, BandWidth, DistN, FragStep);
	Lib.LoadChains(ChainsFN);
	Lib.SetTrainingFrags();
	Lib.AssertSames();
	Lib.Train();
	uint MinClusterSize = Lib.GetMinClusterSize();
	if (MinClusterSize < MinSize)
		Die("MinSize %u", MinClusterSize);
	Lib.AssignReps();
	Lib.LogClusters();
	ProgressLog("MinSize %u (%.2f%%)\n",
		MinSize, GetPct(MinSize, Lib.GetFragCount()));
	Lib.ToSpec(opt(output));
	Lib.AssignLetters();
	Lib.ToFasta(opt(fasta));
	}
