#include "myutils.h"
#include "sort.h"
#include "dbsearcher.h"
#include "chainreader2.h"
#include <thread>

uint GetUBits(const vector<uint> &KmerBitsQ, const vector<uint> &KmerBitsR);

void DBSearcher::USort(const vector<uint> &QueryKmerBits,
  vector<uint> &DBChainIndexes, vector<uint> &Order)
	{
	StartTimer(USort);
	DBChainIndexes.clear();

	const uint MinU = uint(round(m_Params->m_MinU));
	vector<uint> Us;
	uint DBSize = GetDBSize();
	for (uint Idx = 0; Idx < DBSize; ++Idx)
		{
		const vector<uint> &DBKmerBits = *m_DBKmerBitsVec[Idx];
		uint U = GetUBits(QueryKmerBits, DBKmerBits);
		if (U < MinU)
			continue;
		Us.push_back(U);
		DBChainIndexes.push_back(Idx);
		}
	EndTimer(USort);
	const uint N = SIZE(Us);
	if (N > 0)
		{
		Order.resize(N);
		QuickSortOrderDesc(Us.data(), N, Order.data());
		}
	else
		Order.clear();
	}

void DBSearcher::StaticThreadUSort(uint ThreadIndex,
  DBSearcher *ptrDBS, ChainReader2 *ptrQCR)
	{
	ptrDBS->ThreadUSort(ThreadIndex, *ptrQCR);
	}

void DBSearcher::ThreadUSort(uint ThreadIndex, ChainReader2 &QCR)
	{
	asserta(ThreadIndex < SIZE(m_DAs));
	DSSAligner &DA = *m_DAs[ThreadIndex];
	DA.m_Params = m_Params;
	DSS D;
	D.SetParams(*m_Params);
	vector<uint> Order;
	vector<byte> MuLettersQ;
	vector<uint> KmersQ;
	vector<uint> KmerBitsQ;
	vector<vector<byte> > ProfileQ;
	for (;;)
		{
		m_Lock.lock();
		const PDBChain *ptrChainQ = QCR.GetNext();
		if (ptrChainQ == 0)
			{
			m_Lock.unlock();
			break;
			}
		m_Lock.unlock();
		const PDBChain &ChainQ = *ptrChainQ;
		D.Init(ChainQ);
		D.GetMuLetters(MuLettersQ);
		D.GetMuKmers(MuLettersQ, KmersQ);
		D.GetMuKmerBits(KmersQ, KmerBitsQ);
		vector<uint> DBIdxs;
		USort(KmerBitsQ, DBIdxs, Order);
		const uint N = SIZE(Order);
		if (N == 0)
			continue;
		const uint MAXACCEPTS = m_Params->m_MaxAccepts;
		const uint MAXREJECTS = m_Params->m_MaxRejects;
		uint AcceptCount = 0;
		uint RejectCount = 0;
		D.GetProfile(ProfileQ);
		for (uint i = 0; i < N; ++i)
			{
			if (AcceptCount >= MAXACCEPTS)
				break;
			if (RejectCount >= MAXREJECTS)
				break;
			uint ChainIndexR = DBIdxs[Order[i]];
			const PDBChain &ChainR = *m_DBChains[ChainIndexR];
			if (ChainQ.m_Label == ChainR.m_Label)
				continue;
			const vector<byte> &MuLettersR = *m_DBMuLettersVec[ChainIndexR];
			const vector<vector<byte> > &ProfileR = *m_DBProfiles[ChainIndexR];
			DA.Align_MuFilter(ChainQ, ChainR, 
			  MuLettersQ, MuLettersR, ProfileQ, ProfileR);
			if (DA.m_Path.empty())
				continue;
			if (DA.m_EvalueA > m_MaxEvalue)
				{
				++RejectCount;
				continue;
				}
			++AcceptCount;
			if (DA.GetEvalue(true) <= m_MaxEvalue)
				{
				DA.ToTsv(m_fTsv, true);
				BaseOnAln(DA, true);
				}
			}
		}
	}

void DBSearcher::RunUSort(ChainReader2 &QCR)
	{
	m_AlnsPerThreadPerSec = FLT_MAX;
	time_t t_start = time(0);
	m_Secs = UINT_MAX;
	uint ThreadCount = GetRequestedThreadCount();
#if TIMING
	asserta(ThreadCount == 1);
#endif
	m_PairIndex = UINT_MAX;
	m_PairCount = UINT_MAX;
	m_NextChainIndex1 = UINT_MAX;
	m_NextChainIndex2 = UINT_MAX;
	m_NextQueryIdx = UINT_MAX;
	m_NextDBIdx = UINT_MAX;
	m_ProcessedQueryCount = 0;

	m_NextQueryIdx = 0;

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadUSort, ThreadIndex, this, &QCR);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	time_t t_end = time(0);
	m_Secs = uint(t_end - t_start);
	if (m_Secs == 0)
		m_Secs = 1;
	m_AlnsPerThreadPerSec = float(DSSAligner::m_AlnCount)/(m_Secs*ThreadCount);
	}
