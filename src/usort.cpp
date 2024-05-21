#include "myutils.h"
#include "sort.h"
#include "dbsearcher.h"
#include <thread>

uint GetUBits(const vector<uint> &KmerBitsQ, const vector<uint> &KmerBitsR);

void DBSearcher::USort(uint QueryChainIndex, vector<uint> &DBChainIndexes,
  vector<uint> &Order)
	{
	StartTimer(USort);
	DBChainIndexes.clear();
	asserta(m_QueryChainCount + m_DBChainCount == m_ChainCount);
	asserta(QueryChainIndex < m_QueryChainCount);

	const uint MinU = m_Params->m_MinU;
	const vector<uint> &QueryKmerBits = m_KmerBitsVec[QueryChainIndex];
	vector<uint> Us;
	uint DBSize = GetDBSize();
	for (uint Idx = 0; Idx < DBSize; ++Idx)
		{
		uint DBChainIndex = GetDBChainIndex(Idx);
		const vector<uint> &DBKmerBits = m_KmerBitsVec[DBChainIndex];
		uint U = GetUBits(QueryKmerBits, DBKmerBits);
		if (U < MinU)
			continue;
		Us.push_back(U);
		DBChainIndexes.push_back(DBChainIndex);
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

void DBSearcher::StaticThreadUSort(uint ThreadIndex, DBSearcher *ptrDBS)
	{
	ptrDBS->ThreadUSort(ThreadIndex);
	}

void DBSearcher::ThreadUSort(uint ThreadIndex)
	{
	asserta(ThreadIndex < SIZE(m_DAs));
	DSSAligner &DA = *m_DAs[ThreadIndex];
	DA.m_Params = m_Params;
	vector<uint> Order;
	for (;;)
		{
		m_Lock.lock();
		uint ChainIndexQ = m_NextQueryIdx;
		if (ChainIndexQ < m_QueryChainCount)
			{
			ProgressStep(ChainIndexQ, m_QueryChainCount, "Searching (u-sort)");
			++m_NextQueryIdx;
			}
		else if (ChainIndexQ == m_QueryChainCount)
			{
			m_Lock.unlock();
			break;
			}
		else
			asserta(false);
		m_Lock.unlock();
		vector<uint> DBIdxs;
		USort(ChainIndexQ, DBIdxs, Order);
		const uint N = SIZE(Order);
		if (N == 0)
			continue;
		const uint MAXACCEPTS = m_Params->m_MaxAccepts;
		const uint MAXREJECTS = m_Params->m_MaxRejects;
		uint AcceptCount = 0;
		uint RejectCount = 0;
		const PDBChain &ChainQ = *m_Chains[ChainIndexQ];
		const vector<byte> &ComboLettersQ = m_ComboLettersVec[ChainIndexQ];
		const vector<vector<byte> > &ProfileQ = m_Profiles[ChainIndexQ];
		for (uint i = 0; i < N; ++i)
			{
			if (AcceptCount >= MAXACCEPTS)
				break;
			if (RejectCount >= MAXREJECTS)
				break;
			uint ChainIndexR = GetDBChainIndex(DBIdxs[Order[i]]);
			if (ChainIndexQ == ChainIndexR)
				continue;
			const PDBChain &ChainR = *m_Chains[ChainIndexR];
			const vector<byte> &ComboLettersR = m_ComboLettersVec[ChainIndexR];
			const vector<vector<byte> > &ProfileR = m_Profiles[ChainIndexR];
			DA.Align_ComboFilter(ChainQ, ChainR, 
			  ComboLettersQ, ComboLettersR, ProfileQ, ProfileR);
			if (DA.m_PathAB.empty())
				continue;
			if (DA.m_EvalueAB > m_MaxEvalue)
				{
				++RejectCount;
				continue;
				}
			++AcceptCount;
			DA.ToTsv(m_fTsv, m_MaxEvalue);
			OnAln(ChainIndexQ, ChainIndexR, DA);
			}
		}
	}

void DBSearcher::RunUSort()
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
	m_ProcessedQueryCount = UINT_MAX;

	m_NextQueryIdx = 0;

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadUSort, ThreadIndex, this);
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
