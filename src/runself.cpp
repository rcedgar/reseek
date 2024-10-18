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
	const bool HasSelfRevScores = !m_DBSelfRevScores.empty();
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
			const vector<vector<byte> > *ptrProfile1 = m_DBProfiles[ChainIndex1];
			const vector<byte> *ptrComboLetters1 = (m_DBComboLettersVec.empty() ? 0 : m_DBComboLettersVec[ChainIndex1]);
			const vector<uint> *ptrKmerBits1 = (m_DBKmerBitsVec.empty() ? 0 : m_DBKmerBitsVec[ChainIndex1]);
			float SelfRevScore1 = HasSelfRevScores ? m_DBSelfRevScores[ChainIndex1] : FLT_MAX;
			DA.SetQuery(Chain1, ptrProfile1, ptrKmerBits1, ptrComboLetters1, SelfRevScore1);
			}

		const PDBChain &Chain2 = *m_DBChains[ChainIndex2];
		const vector<vector<byte> > *ptrProfile2 = m_DBProfiles[ChainIndex2];
		const vector<byte> *ptrComboLetters2 = (m_DBComboLettersVec.empty() ? 0 : m_DBComboLettersVec[ChainIndex2]);
		const vector<uint> *ptrKmerBits2 = (m_DBKmerBitsVec.empty() ? 0 : m_DBKmerBitsVec[ChainIndex2]);
		float SelfRevScore2 = HasSelfRevScores ? m_DBSelfRevScores[ChainIndex2] : FLT_MAX;
		DA.SetTarget(Chain2, ptrProfile2, ptrKmerBits2, ptrComboLetters2, SelfRevScore2);

		DA.AlignQueryTarget();
		if (!DA.m_Path.empty())
			{
			BaseOnAln(DA, true);
			if (ChainIndex1 != ChainIndex2)
				BaseOnAln(DA, false);
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
		++m_ProcessedQueryCount;
		++m_NextChainIndex1;
		m_NextChainIndex2 = m_NextChainIndex1;
		}
	++m_PairIndex;
	++m_ProcessedPairCount;
	m_Lock.unlock();
	return true;
	}

void DBSearcher::RunSelf()
	{
	m_D.SetParams(*m_Params);
	for (uint i = 0; i < SIZE(m_DAs); ++i)
		m_DAs[i]->m_Params = m_Params;

	asserta(!m_Params->m_USort);

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

	uint ChainCount = GetDBChainCount();
	m_PairIndex = 0;
	m_PairCount = ChainCount + (ChainCount*(ChainCount-1))/2;
	m_NextChainIndex1 = 0;
	m_NextChainIndex2 = 0;

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
