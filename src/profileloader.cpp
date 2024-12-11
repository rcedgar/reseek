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

	const uint L = Chain.GetSeqLength();
	if (L > 512)
		return 0;

	PDBChain RevChain = Chain;
	RevChain.Reverse();
	vector<vector<byte> > RevProfile;
	D.Init(RevChain);
	D.GetProfile(RevProfile);
	DA.SetQuery(Chain, &Profile, 0, 0, 0, 0);
	DA.SetTarget(RevChain, &RevProfile, 0, 0, 0, 0);
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

	DA.SetParams(DA_Params);

	for (;;)
		{
		PDBChain *Chain = m_CR->GetNext();
		if (Chain == 0)
			return;
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
		vector<byte> *MuLetters = m_MuLetters == 0 ? 0 : new vector<byte>;
		vector<uint> *Kmers = m_MuLetters == 0 ? 0 : new vector<uint>;
		float SelfRevScore = FLT_MAX;

		D.Init(*Chain);
		if (m_Profiles != 0) D.GetProfile(*ptrProfile);
		if (m_MuLetters != 0) D.GetMuLetters(*MuLetters);
		if (m_KmerBitsVec != 0) D.GetMuKmerBits(*Kmers, *KmerBits);
		if (m_MuLetters != 0) D.GetMuKmers(*MuLetters, *Kmers);
		if (m_SelfRevScores != 0) SelfRevScore = GetSelfRevScore(DA, D, *Chain, *ptrProfile);

		m_Lock.lock();
		Chain->m_Idx = SIZE(*m_Chains);
		if (m_Chains != 0) m_Chains->push_back(Chain);
		if (m_Profiles != 0) m_Profiles->push_back(ptrProfile);
		if (m_KmerBitsVec != 0) m_KmerBitsVec->push_back(KmerBits);
		if (m_MuLetters != 0) m_MuLetters->push_back(MuLetters);
		if (m_KmersVec != 0) m_KmersVec->push_back(Kmers);
		if (m_SelfRevScores != 0) m_SelfRevScores->push_back(SelfRevScore);
		m_Lock.unlock();
		}
	}

void ProfileLoader::Load(
  const DSSParams &Params,
  ChainReader2 &CR,
  vector<PDBChain *> *Chains,
  vector<vector<vector<byte> > *> *Profiles,
  vector<vector<byte> *> *MuLetters,
  vector<vector<uint> *> *KmersVec,
  vector<vector<uint> *> *KmerBitsVec,
  vector<float> *SelfRevScores,
  uint ThreadCount)
	{
	m_MinChainLength = 1;
	if (optset_minchainlength)
		m_MinChainLength = opt_minchainlength;

	if (KmerBitsVec != 0)
		asserta(MuLetters != 0);

	m_Params = &Params;
	m_CR = &CR;
	m_Chains = Chains;
	m_Profiles = Profiles;
	m_KmerBitsVec = KmerBitsVec;
	m_KmersVec = KmersVec;
	m_MuLetters = MuLetters;
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
