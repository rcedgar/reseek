#include "myutils.h"
#include "dbsearcher.h"
#include "chainreader2.h"
#include "timing.h"

void DBSearcher::StaticThreadBodyQuery(uint ThreadIndex, DBSearcher *ptrDBS,
  ChainReader2 *ptrQueryCR)
	{
	ptrDBS->ThreadBodyQuery(ThreadIndex, ptrQueryCR);
	}

void DBSearcher::ThreadBodyQuery(uint ThreadIndex, ChainReader2 *ptrQueryCR)
	{
	asserta(ThreadIndex < SIZE(m_DAs));
	DSSAligner &DA = *m_DAs[ThreadIndex];
	const uint DBChainCount = GetDBChainCount();

	vector<vector<byte> > Profile1;
	vector<byte> ComboLetters1;
	vector<uint> Kmers1;
	vector<uint> KmerBits1;

	DSS D;
	D.SetParams(*m_Params);

	for (;;)
		{
		PDBChain *Chain1 = ptrQueryCR->GetNext();
		if (Chain1 == 0)
			return;
		D.Init(*Chain1);
		D.GetProfile(Profile1);
		if (m_Params->m_Omega > 0)
			D.GetComboLetters(ComboLetters1);
		if (m_Params->m_MinU > 0)
			{
			D.GetComboKmers(ComboLetters1, Kmers1);
			D.GetComboKmerBits(Kmers1, KmerBits1);
			}

		const vector<byte> *ptrComboLetters1 = (ComboLetters1.empty() ? 0 : &ComboLetters1);
		const vector<uint> *ptrKmerBits1 = (KmerBits1.empty() ? 0 : &KmerBits1);
		DA.SetQuery(*Chain1, &Profile1, ptrKmerBits1, ptrComboLetters1);

		for (uint DBChainIdx = 0; DBChainIdx < DBChainCount; ++DBChainIdx)
			{
			const PDBChain &Chain2 = *m_DBChains[DBChainIdx];
			const vector<vector<byte> > *ptrProfile2 = m_DBProfiles[DBChainIdx];
			const vector<byte> *ptrComboLetters2 = (m_DBComboLettersVec.empty() ? 0 : m_DBComboLettersVec[DBChainIdx]);
			const vector<uint> *ptrKmerBits2 = (m_DBKmerBitsVec.empty() ? 0 : m_DBKmerBitsVec[DBChainIdx]);
			DA.SetTarget(Chain2, ptrProfile2, ptrKmerBits2, ptrComboLetters2);

			DA.AlignQueryTarget();
			if (!DA.m_Path.empty())
				BaseOnAln(DA, true);
			++m_ProcessedPairCount;
			}
		delete Chain1;
		m_Lock.lock();
		++m_ProcessedQueryCount;
		time_t Now = time(0);
		if (Now > m_LastProgress)
			{
			Progress("%s query chains\r", IntToStr(m_ProcessedQueryCount));
			m_LastProgress = Now;
			}
		m_Lock.unlock();
		}
	}

void DBSearcher::RunQuery(ChainReader2 &QCR)
	{
	m_D.SetParams(*m_Params);
	for (uint i = 0; i < SIZE(m_DAs); ++i)
		m_DAs[i]->m_Params = m_Params;

	if (m_Params->m_USort)
		{
		RunUSort(QCR);
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
	m_ProcessedQueryCount = 0;

	uint ChainCount = GetDBChainCount();
	m_PairIndex = 0;
	m_PairCount = (ChainCount*(ChainCount-1))/2;
	m_NextChainIndex1 = 0;
	m_NextChainIndex2 = 1;

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadBodyQuery, ThreadIndex, this, &QCR);
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
	Progress("%s query chains\n", IntToStr(m_ProcessedQueryCount));
	RunStats();
	}
