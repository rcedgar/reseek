#include "myutils.h"
#include "dbsearcher.h"
#include "chainreader2.h"
#include "timing.h"
#include "mx.h"

float GetSelfRevScore(DSSAligner &DA, const PDBChain &Chain,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers,
					  const vector<vector<byte> > &Profile);

void DBSearcher::StaticThreadBodyQuery(uint ThreadIndex, DBSearcher *ptrDBS,
  ChainReader2 *ptrQueryCR)
	{
	ptrDBS->ThreadBodyQuery(ThreadIndex, ptrQueryCR);
	}

//float DBSearcher::GetSelfRevScore(const PDBChain &Chain,
//  const vector<vector<byte> > &Profile, DSSAligner &DA)
//	{
//	if (opt_selfrev0)
//		return 0;
//
//	DA.SetQuery(Chain, &Profile, 0, 0, FLT_MAX);
//	DA.Align_QRev();
//	return DA.m_AlnFwdScore;
//	}

void DBSearcher::ThreadBodyQuery(uint ThreadIndex, ChainReader2 *ptrQueryCR)
	{
	asserta(ThreadIndex < SIZE(m_DAs));
	DSSAligner &DA = *m_DAs[ThreadIndex];
	const uint DBChainCount = GetDBChainCount();

	vector<vector<byte> > Profile1;
	vector<byte> MuLetters1;
	vector<uint> MuKmers1;

	DSS D;
	D.SetParams(*m_Params);

	for (;;)
		{
		PDBChain *Chain1 = ptrQueryCR->GetNext();
		if (Chain1 == 0)
			{
			static bool ExitMsgDone = false;
			m_Lock.lock();
			if (!ExitMsgDone)
				{
				ProgressLog("Before exit threads %.4g\n", GetMemUseBytes());
				ExitMsgDone = true;
				}
			m_Lock.unlock();
			return;
			}
		D.Init(*Chain1);
		D.GetProfile(Profile1);
		D.GetMuLetters(MuLetters1);
		D.GetMuKmers(MuLetters1, MuKmers1);

		const vector<byte> *ptrMuLetters1 = (MuLetters1.empty() ? 0 : &MuLetters1);
		const vector<uint> *ptrMuKmers1 = (MuKmers1.empty() ? 0 : &MuKmers1);
		float SelfRevScore =
			GetSelfRevScore(DA, *Chain1, ptrMuLetters1, ptrMuKmers1, Profile1);
		DA.SetQuery(*Chain1, &Profile1, ptrMuLetters1, ptrMuKmers1, SelfRevScore);

		for (uint DBChainIdx = 0; DBChainIdx < DBChainCount; ++DBChainIdx)
			{
			m_Lock.lock();

			if (DBChainIdx%10000 == 0)
				{
				time_t now = time(0);
				if (now > m_LastProgress)
					{
					Progress("%s query chains\r",
							 IntToStr(m_ProcessedQueryCount));
					m_LastProgress = now;
					}
				}
			m_Lock.unlock();
			const PDBChain &Chain2 = *m_DBChains[DBChainIdx];
			if (opt_noself && Chain1->m_Label == Chain2.m_Label)
				continue;
			const vector<vector<byte> > *ptrProfile2 = m_DBProfiles[DBChainIdx];
			const vector<byte> *ptrMuLetters2 = (m_DBMuLettersVec.empty() ? 0 : m_DBMuLettersVec[DBChainIdx]);
			const vector<uint> *ptrMuKmers2 = (m_DBMuLettersVec.empty() ? 0 : m_DBMuKmersVec[DBChainIdx]);
			float SelfRevScore = (m_DBSelfRevScores.empty() ? FLT_MAX : m_DBSelfRevScores[DBChainIdx]);
			DA.SetTarget(Chain2, ptrProfile2, ptrMuLetters2, ptrMuKmers2, SelfRevScore);
			DA.AlignQueryTarget();
			if (!DA.m_Path.empty())
				BaseOnAln(DA, true);
			++m_ProcessedPairCount;
			}
		DA.UnsetQuery();
		delete Chain1;
		++m_ProcessedQueryCount;
		}
	}

void DBSearcher::RunQuery(ChainReader2 &QCR)
	{
	for (uint i = 0; i < SIZE(m_DAs); ++i)
		m_DAs[i]->SetParams(*m_Params);

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
	ProgressLog("After exit threads %.4g\n", GetMemUseBytes());

	time_t t_end = time(0);
	m_Secs = uint(t_end - t_start);
	if (m_Secs == 0)
		m_Secs = 1;
	m_AlnsPerThreadPerSec = float(DSSAligner::m_AlnCount)/(m_Secs*ThreadCount);
	Progress("%s query chains\n", IntToStr(m_ProcessedQueryCount));
	RunStats();
	}
