#include "myutils.h"
#include "chainreader2.h"
#include "dss.h"
#include "dssaligner.h"
#include "profileloader.h"

float GetSelfRevScore(DSSAligner &DA, const PDBChain &Chain,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers,
					  const vector<vector<byte> > &Profile)
	{
	if (opt_selfrev0)
		return 0;

	DA.SetQuery(Chain, &Profile, ptrMuLetters, ptrMuKmers, FLT_MAX);

	PDBChain ChainRev;
	Chain.GetReverse(ChainRev);
	vector<vector<byte> > ProfileRev;

	const uint L = Chain.GetSeqLength();
	const uint FeatureCount = SIZE(Profile);
	ProfileRev.resize(FeatureCount);
	for (uint F = 0; F < FeatureCount; ++F)
		{
		const vector<byte> &Row = Profile[F];
		vector<byte> &RowRev = ProfileRev[F];
		RowRev.reserve(L);
		for (uint i = 0; i < L; ++i)
			RowRev.push_back(Row[L-i-1]);
		}

#if DEBUG
	{
	const uint L = Chain.GetSeqLength();
	asserta(L > 0);
	asserta(SIZE(ProfileRev) == FeatureCount);
	for (uint F = 0; F < FeatureCount; ++F)
		{
		const vector<byte> &Row = Profile[F];
		const vector<byte> &RowRev = ProfileRev[F];
		asserta(SIZE(Row) == L);
		asserta(SIZE(RowRev) == L);
		for (uint Pos = 0; Pos < L; ++Pos)
			asserta(Row[Pos] == RowRev[L-Pos-1]);
		}
	}
#endif

	vector<byte> MuLettersRev;
	vector<uint> MuKmersRev;

	const vector<byte> *ptrMuLettersRev = 0;
	const vector<uint> *ptrMuKmersRev = 0;

	if (ptrMuLetters != 0)
		{
		asserta(ptrMuKmers != 0);
		const vector<byte> &MuLetters = *ptrMuLetters;
		const vector<uint> &MuKmers = *ptrMuKmers;

		const uint L = SIZE(MuLetters);
		MuLettersRev.reserve(L);
		for (uint i = 0; i < L; ++i)
			MuLettersRev.push_back(MuLetters[L-i-1]);

		const uint K = SIZE(MuKmers);
		MuKmersRev.reserve(K);
		for (uint i = 0; i < K; ++i)
			MuKmersRev.push_back(MuKmers[K-i-1]);


		ptrMuLettersRev = &MuLettersRev;
		ptrMuKmersRev = &MuKmersRev;
		}
	DA.SetTarget(ChainRev, &ProfileRev, ptrMuLettersRev, ptrMuKmersRev, FLT_MAX);
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
	DA_Params.m_OwnScoreMxs = false;

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
		vector<byte> *MuLetters = m_MuLetters == 0 ? 0 : new vector<byte>;
		vector<uint> *MuKmers = m_MuLetters == 0 ? 0 : new vector<uint>;
		float SelfRevScore = FLT_MAX;

		D.Init(*Chain);
		if (m_Profiles != 0) D.GetProfile(*ptrProfile);
		if (m_MuLetters != 0) D.GetMuLetters(*MuLetters);
		if (m_MuLetters != 0) D.GetMuKmers(*MuLetters, *MuKmers);
		if (m_SelfRevScores != 0) SelfRevScore =
			GetSelfRevScore(DA, *Chain, MuLetters, MuKmers, *ptrProfile);

		m_Lock.lock();
		Chain->m_Idx = SIZE(*m_Chains);
		if (m_Chains != 0) m_Chains->push_back(Chain);
		if (m_Profiles != 0) m_Profiles->push_back(ptrProfile);
		if (m_MuLetters != 0) m_MuLetters->push_back(MuLetters);
		if (m_KmersVec != 0) m_KmersVec->push_back(MuKmers);
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
  vector<float> *SelfRevScores,
  uint ThreadCount)
	{
	m_MinChainLength = 1;
	if (optset_minchainlength)
		m_MinChainLength = opt_minchainlength;

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

	Progress("%s target chains loaded\n", IntToStr(m_Count));
	}
