#include "myutils.h"
#include "dbsearcher.h"
#include "timing.h"

void DBSearcher::StaticThreadBodySelf(uint ThreadIndex, DBSearcher *ptrDBS)
	{
	ptrDBS->ThreadBodySelf(ThreadIndex);
	}

void DBSearcher::ThreadBodySelf(uint ThreadIndex)
	{
	asserta(ThreadIndex < SIZE(m_DAs));
	uint PrevChainIndex1 = UINT_MAX;
	DSSAligner &DA = *m_DAs[ThreadIndex];
	for (;;)
		{
		uint ChainIndex1, ChainIndex2;
		bool Ok = GetNextPairSelf(ChainIndex1, ChainIndex2);
		if (!Ok)
			break;

		if (ChainIndex1 == PrevChainIndex1)
			++m_QPCacheHits;
		else
			{
			++m_QPCacheMisses;
			const PDBChain &Chain1 = *m_DBChains[ChainIndex1];
			const vector<vector<byte> > &Profile1 = *m_DBProfiles[ChainIndex1];
			const vector<byte> &ComboLetters1 = *m_DBComboLettersVec[ChainIndex1];
			const vector<uint> &KmerBits1 = *m_DBKmerBitsVec[ChainIndex1];
			DA.SetQuery(Chain1, &Profile1, &KmerBits1, &ComboLetters1);
			}

		const PDBChain &Chain2 = *m_DBChains[ChainIndex2];
		const vector<vector<byte> > &Profile2 = *m_DBProfiles[ChainIndex2];
		const vector<byte> &ComboLetters2 = *m_DBComboLettersVec[ChainIndex2];
		const vector<uint> &KmerBits2 = *m_DBKmerBitsVec[ChainIndex2];
		DA.SetTarget(Chain2, &Profile2, &KmerBits2, &ComboLetters2);

		DA.AlignQueryTarget();
		if (!DA.m_Path.empty())
			{
			if (m_CollectTestStats)
				{
				m_Lock.lock();
				asserta(ChainIndex1 < SIZE(m_TestStatsVec));
				vector<float> &v1 = m_TestStatsVec[ChainIndex1];
				v1.push_back(DA.m_TestStatisticA);
				m_Lock.unlock();
				}

			BaseOnAln(ChainIndex1, ChainIndex2, DA, true);
			if (m_QuerySelf)
				{
				if (m_CollectTestStats)
					{
					m_Lock.lock();
					asserta(ChainIndex2 < SIZE(m_TestStatsVec));
					vector<float> &v2 = m_TestStatsVec[ChainIndex2];
					v2.push_back(DA.m_TestStatisticB);
					m_Lock.unlock();
					}

				BaseOnAln(ChainIndex1, ChainIndex2, DA, false);
				}
			}
		PrevChainIndex1 = ChainIndex1;
		}
	}

bool DBSearcher::GetNextPairSelf(uint &ChainIndex1, uint &ChainIndex2)
	{
	ChainIndex1 = UINT_MAX;
	ChainIndex2 = UINT_MAX;
	m_Lock.lock();
	if (m_PairIndex == m_PairCount)
		{
		m_Lock.unlock();
		return false;
		}
	ProgressStep(m_PairIndex, m_PairCount, "Aligning");
	ChainIndex1 = m_NextChainIndex1;
	ChainIndex2 = m_NextChainIndex2;
	uint ChainCount = GetDBChainCount();
	assert(m_NextChainIndex1 < ChainCount);
	assert(m_NextChainIndex2 < ChainCount);
	++m_NextChainIndex2;
	if (m_NextChainIndex2 == ChainCount)
		{
		++m_NextChainIndex1;
		m_NextChainIndex2 = m_NextChainIndex1 + 1;
		}
	++m_PairIndex;
	m_Lock.unlock();
	return true;
	}

void DBSearcher::RunSelf()
	{
	m_D.SetParams(*m_Params);
	for (uint i = 0; i < SIZE(m_DAs); ++i)
		m_DAs[i]->m_Params = m_Params;

	if (m_Params->m_USort)
		{
		RunUSort();
		return;
		}

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

	uint ChainCount = GetDBChainCount();
	m_PairIndex = 0;
	m_PairCount = (ChainCount*(ChainCount-1))/2;
	m_NextChainIndex1 = 0;
	m_NextChainIndex2 = 1;

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadBodySelf, ThreadIndex, this);
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
	RunStats();
	}
