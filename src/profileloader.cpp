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
	D.SetParams(*m_Params);

	vector<uint> Kmers;
	for (;;)
		{
		PDBChain *Chain = m_CR->GetNext();
		if (Chain == 0)
			return;
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

		D.Init(*Chain);
		if (m_Profiles != 0) D.GetProfile(*ptrProfile);
		if (m_ComboLetters != 0) D.GetComboLetters(*ComboLetters);
		if (m_ComboLetters != 0) D.GetComboKmers(*ComboLetters, Kmers);
		if (m_KmerBitsVec != 0) D.GetComboKmerBits(Kmers, *KmerBits);

		m_Lock.lock();
		if (m_Chains != 0) m_Chains->push_back(Chain);
		if (m_Profiles != 0) m_Profiles->push_back(ptrProfile);
		if (m_KmerBitsVec != 0) m_KmerBitsVec->push_back(KmerBits);
		if (m_ComboLetters != 0) m_ComboLetters->push_back(ComboLetters);
		m_Lock.unlock();
		}
	}

void ProfileLoader::Load(
  ChainReader2 &CR,
  vector<PDBChain *> *Chains,
  vector<vector<vector<byte> > *> *Profiles,
  vector<vector<byte> *> *ComboLetters,
  vector<vector<uint> *> *KmerBitsVec,
  const DSSParams &Params,
  uint ThreadCount)
	{
	if (KmerBitsVec != 0)
		asserta(ComboLetters != 0);

	m_CR = &CR;
	m_Chains = Chains;
	m_Profiles = Profiles;
	m_KmerBitsVec = KmerBitsVec;
	m_ComboLetters = ComboLetters;
	m_Params = &Params;

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

void cmd_test()
	{
	ProfileLoader PL;
	ChainReader2 CR;
	CR.Open(g_Arg1);

	DSSParams Params;
	Params.SetFromCmdLine(10000);

	vector<PDBChain *> Chains;
	vector<vector<vector<byte> > *> Profiles;
	vector<vector<byte> *> ComboLetters;
	vector<vector<uint> *> KmersVec;
	vector<vector<uint> *> KmerBitsVec;

	uint ThreadCount = GetRequestedThreadCount();
	PL.Load(CR, &Chains, &Profiles, &ComboLetters,
	  &KmerBitsVec, Params, ThreadCount);
	ProgressLog("Done.\n");
	}
