#include "myutils.h"
#include "chainreader2.h"
#include "dss.h"
#include "dssaligner.h"
#include "profileloader.h"

float GetSelfRevScore(DSSAligner &DA, DSS &D, 
  const PDBChain &Chain, const vector<vector<byte> > &Profile)
	{
	if (opt_selfrev0)
		return 0;

	PDBChain RevChain = Chain;
	RevChain.Reverse();
	vector<vector<byte> > RevProfile;
	D.Init(RevChain);
	D.GetProfile(RevProfile);
	DA.SetQuery(Chain, &Profile, 0, 0, 0);
	DA.SetTarget(RevChain, &RevProfile, 0, 0, 0);
	DA.AlignQueryTarget();
	return DA.m_AlnFwdScore;
	}

void ProfileLoader::StaticThreadBody(uint ThreadIndex, ProfileLoader *PL)
	{
	PL->ThreadBody(ThreadIndex);
	}

void ProfileLoader::ThreadBody(uint ThreadIndex)
	{
	DSS D;
	D.SetParams(*m_Params);

	DSSAligner DA;
	DSSParams DA_Params = *m_Params;
	DA_Params.m_UsePara = false;
	DA_Params.m_Omega = 0;

	DA.m_Params = &DA_Params;

	vector<uint> Kmers;
	for (;;)
		{
		PDBChain *Chain = m_CR->GetNext();
		if (Chain == 0)
			{
			//static mutex m;
			//m.lock();
			//ProgressLogPrefix("Chain=0 %.3g\n", GetMemUseBytes());
			//m.unlock();
			return;
			}
		if (Chain->GetSeqLength() < m_MinChainLength)
			{
			delete Chain;
			continue;
			}
		time_t Now = time(0);
		m_Lock.lock();
		++m_Count;
		if (Now > m_LastProgress)
			{
			Progress("%s target chains loaded\r", IntToStr(m_Count));
			m_LastProgress = Now;
			}
		m_Lock.unlock();

		vector<vector<byte> > *ptrProfile = m_Profiles == 0 ? 0 : new vector<vector<byte> >;
		vector<uint> *KmerBits = m_KmerBitsVec == 0 ? 0 : new vector<uint>;
		vector<byte> *ComboLetters = m_ComboLetters == 0 ? 0 : new vector<byte>;
		float SelfRevScore = FLT_MAX;

		D.Init(*Chain);
		if (m_Profiles != 0) D.GetProfile(*ptrProfile);
		if (m_ComboLetters != 0) D.GetComboLetters(*ComboLetters);
		if (m_ComboLetters != 0) D.GetComboKmers(*ComboLetters, Kmers);
		if (m_KmerBitsVec != 0) D.GetComboKmerBits(Kmers, *KmerBits);
		if (m_SelfRevScores != 0) SelfRevScore = GetSelfRevScore(DA, D, *Chain, *ptrProfile);

		m_Lock.lock();
		Chain->m_Idx = SIZE(*m_Chains);
		if (m_Chains != 0) m_Chains->push_back(Chain);
		if (m_Profiles != 0) m_Profiles->push_back(ptrProfile);
		if (m_KmerBitsVec != 0) m_KmerBitsVec->push_back(KmerBits);
		if (m_ComboLetters != 0) m_ComboLetters->push_back(ComboLetters);
		if (m_SelfRevScores != 0) m_SelfRevScores->push_back(SelfRevScore);
		m_Lock.unlock();
		}
	}

void ProfileLoader::Load(
  const DSSParams &Params,
  ChainReader2 &CR,
  vector<PDBChain *> *Chains,
  vector<vector<vector<byte> > *> *Profiles,
  vector<vector<byte> *> *ComboLetters,
  vector<vector<uint> *> *KmerBitsVec,
  vector<float> *SelfRevScores,
  uint ThreadCount)
	{
	m_MinChainLength = 1;
	if (optset_minchainlength)
		m_MinChainLength = opt_minchainlength;

	if (KmerBitsVec != 0)
		asserta(ComboLetters != 0);

	m_Params = &Params;
	m_CR = &CR;
	m_Chains = Chains;
	m_Profiles = Profiles;
	m_KmerBitsVec = KmerBitsVec;
	m_ComboLetters = ComboLetters;
	m_SelfRevScores = SelfRevScores;

	Chains->clear();
	Profiles->clear();

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadBody, ThreadIndex, this);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();

	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	Progress("%s target chains loaded\n", IntToStr(m_Count));
	}
