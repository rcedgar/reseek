#include "myutils.h"
#include "dbsearcher.h"

float GetSelfRevScore(DSSAligner &DA, DSS &D, const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers);

static DBSearcher *s_DBS;
static uint s_ChainCount;
static uint s_NextChainIdx;
static mutex s_Lock;

static void ThreadBody(uint ThreadIndex)
	{
	DSS D;
	DSSAligner DA;

	for (;;)
		{
		uint MyChainIdx = UINT_MAX;
		s_Lock.lock();
		if (s_NextChainIdx < s_ChainCount)
			MyChainIdx = s_NextChainIdx++;
		s_Lock.unlock();
		if (MyChainIdx == UINT_MAX)
			return;

		PDBChain *Chain = s_DBS->m_DBChains[MyChainIdx];

		vector<vector<byte> > *ptrProfile = s_DBS->m_DBProfiles[MyChainIdx];
		vector<byte> *MuLetters = s_DBS->m_DBMuLettersVec[MyChainIdx];
		vector<uint> *MuKmers = s_DBS->m_DBMuKmersVec[MyChainIdx];
		D.Init(*Chain);
		float SelfRevScore =
			GetSelfRevScore(DA, D, *Chain, *ptrProfile, MuLetters, MuKmers);
		s_DBS->m_DBSelfRevScores[MyChainIdx] = SelfRevScore;
		}
	}

void DBSearcher::SetSelfRevScores()
	{
	const uint ChainCount = GetDBChainCount();
	m_DBSelfRevScores.clear();
	m_DBSelfRevScores.resize(ChainCount, FLT_MAX);
	s_DBS = this;
	s_ChainCount = GetDBChainCount();
	s_NextChainIdx = 0;
	vector<thread *> ts;
	const uint ThreadCount = GetRequestedThreadCount();
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
