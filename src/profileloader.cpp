#include "myutils.h"
#include "chainreader2.h"
#include "dss.h"
#include "dssaligner.h"
#include "profileloader.h"

void ProfileLoader::StaticThreadBody(uint ThreadIndex, ProfileLoader *PL)
	{
	PL->ThreadBody(ThreadIndex);
	}

void ProfileLoader::ThreadBody(uint ThreadIndex)
	{
	DSS D;
	DSSAligner DA;
	DSSParams DA_Params = *m_Params;
	DA_Params.m_OwnScoreMxs = false;
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
		StartTimer(ProfileLoader_Lock);
		m_Lock.lock();
		++m_Count;
		if (Now > m_LastProgress)
			{
			Progress("%s loaded\r", IntToStr(m_Count));
			m_LastProgress = Now;
			}
		m_Lock.unlock();
		EndTimer(ProfileLoader_Lock);

		StartTimer(ProfileLoader_New);
		vector<vector<byte> > *ptrProfile = m_Profiles == 0 ? 0 : new vector<vector<byte> >;
		vector<byte> *MuLetters = m_MuLetters == 0 ? 0 : new vector<byte>;
		vector<uint> *MuKmers = m_MuLetters == 0 ? 0 : new vector<uint>;
		float SelfRevScore = FLT_MAX;
		EndTimer(ProfileLoader_New);

		D.Init(*Chain, *m_Params);
		if (m_Profiles != 0) D.GetProfile(*ptrProfile);
		if (m_MuLetters != 0) D.GetMuLetters(*MuLetters);
		if (m_MuLetters != 0) D.GetMuKmers(*MuLetters, *MuKmers, m_Params->m_MKFPatternStr);
		if (m_SelfRevScores != 0) SelfRevScore =
			GetSelfRevScore(DA, DA_Params, *Chain, *ptrProfile, MuLetters, MuKmers);

		StartTimer(ProfileLoader_Lock);
		m_Lock.lock();
		Chain->m_Idx = SIZE(*m_Chains);
		if (m_Chains != 0) m_Chains->push_back(Chain);
		if (m_Profiles != 0) m_Profiles->push_back(ptrProfile);
		if (m_MuLetters != 0) m_MuLetters->push_back(MuLetters);
		if (m_KmersVec != 0) m_KmersVec->push_back(MuKmers);
		if (m_SelfRevScores != 0) m_SelfRevScores->push_back(SelfRevScore);
		m_Lock.unlock();
		EndTimer(ProfileLoader_Lock);
		}
	}

void ProfileLoader::Load(
  const DSSParams &Params,
  ChainReader2 &CR,
  vector<PDBChain *> *Chains,
  vector<vector<vector<byte> > *> *Profiles,
  vector<vector<byte> *> *MuLetters,
  vector<vector<uint> *> *KmersVec,
  vector<float> *SelfRevScores,
  uint ThreadCount)
	{
	m_MinChainLength = 1;
	if (optset_minchainlength)
		m_MinChainLength = opt(minchainlength);

	m_Params = &Params;
	m_CR = &CR;
	m_Chains = Chains;
	m_Profiles = Profiles;
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

	Progress("%s chains loaded\n", IntToStr(m_Count));
	}
