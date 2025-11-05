#include "myutils.h"
#include "dbsearcher.h"

float GetSelfRevScore(DSSAligner &DA, DSS &D, const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers)
	{
	if (opt(selfrev0))
		return 0;
	PDBChain RevChain;
	Chain.GetReverse(RevChain);
	vector<vector<byte> > RevProfile;
	D.Init(RevChain);
	D.GetProfile(RevProfile);
	DA.SetQuery(Chain, &Profile, ptrMuLetters, ptrMuKmers, FLT_MAX);
	DA.SetTarget(RevChain, &RevProfile, ptrMuLetters, ptrMuKmers, FLT_MAX);
	DA.AlignQueryTarget();
	return DA.m_AlnFwdScore;
	}

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

#if 0
void cmd_test()
	{
	const string &FN = g_Arg1;

	optset_fast = true;
	opt(fast) = true;
	DSSParams::Init(DM_AlwaysSensitive);

	DBSearcher DBS;
	DBS.LoadDB(FN);
	const uint N = 1000;
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Selfrev");
		DBS.SetSelfRevScores();
		}
	}
#endif