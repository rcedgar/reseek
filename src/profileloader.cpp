#include "myutils.h"
#include "chainreader2.h"
#include "dss.h"
#include "profileloader.h"

void ProfileLoader::StaticThreadBody(uint ThreadIndex, ProfileLoader *PL)
	{
	PL->ThreadBody(ThreadIndex);
	}

void ProfileLoader::ThreadBody(uint ThreadIndex)
	{
	DSS D;
	D.m_Params = m_Params;

	PDBChain *Chain = m_CR->GetNext();
	if (Chain == 0)
		return;
	D.Init(*Chain);
	vector<vector<byte> > *ptrProfile = new vector<vector<byte> >;
	D.GetProfile(*ptrProfile);

	m_Lock.lock();
	m_Chains->push_back(Chain);
	m_Profiles->push_back(ptrProfile);
	m_Lock.unlock();
	}

void ProfileLoader::Load(
  ChainReader2 &CR,
  uint ReserveSize,
  vector<PDBChain *> &Chains,
  vector<vector<vector<byte> > *> &Profiles,
  const DSSParams &Params,
  uint ThreadCount)
	{
	m_CR = &CR;
	m_Chains = &Chains;
	m_Profiles = &Profiles;
	m_Params = &Params;

	Chains.clear();
	Profiles.clear();
	Chains.reserve(ReserveSize);
	Profiles.reserve(ReserveSize);

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
	}
