#include "myutils.h"
#include "cigar.h"
#include "dss.h"
#include "dbsearcher.h"
#include "scop40bench.h"
#include <thread>
#include "timing.h"

void DBSearcher::SetKmersVec()
	{
	StartTimer(SetKmersVec);
	DSS &D = m_D;
	m_KmersVec.clear();
	m_KmerBitsVec.clear();
	asserta(m_ChainCount < 100000);
	m_KmersVec.resize(m_ChainCount);
	m_KmerBitsVec.resize(m_ChainCount);
	for (uint ChainIndex = 0; ChainIndex < m_ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, m_ChainCount, "Set k-mers");
		const PDBChain &Chain = *m_Chains[ChainIndex];
		D.Init(Chain);
		vector<uint> &Kmers = m_KmersVec[ChainIndex];
		vector<uint> &Bits = m_KmerBitsVec[ChainIndex];
		D.GetComboKmers(Kmers);
		D.GetComboKmerBits(Kmers, Bits);
		}
	EndTimer(SetKmersVec);
	}

// Query chains first, then DB chains (unless self)
void DBSearcher::LoadChains(const string &QueryCalFileName, 
  const string &DBCalFileName)
	{
	m_ChainCount = 0;
	m_Chains.clear();
	m_Profiles.clear();

	::ReadChains(QueryCalFileName, m_QueryChains);
	m_QueryChainCount = SIZE(m_QueryChains);

	if (DBCalFileName == "")
		{
		m_QuerySelf = true;
		m_DBChains.clear();
		m_DBChainCount = 0;
		}
	else
		{
		m_QuerySelf = false;
		::ReadChains(DBCalFileName, m_DBChains);
		m_DBChainCount = SIZE(m_DBChains);
		}

	m_ChainCount = m_QueryChainCount + m_DBChainCount;

	m_Chains = m_QueryChains;
	m_Chains.insert(m_Chains.end(),
	  m_DBChains.begin(), m_DBChains.end());
	asserta(SIZE(m_Chains) == m_ChainCount);

	if (m_DBChainCount > 0)
		ProgressLog("%u query chains, %u database chains\n",
		  m_QueryChainCount, m_DBChainCount);
	else
		ProgressLog("All-vs-all %u chains\n", m_QueryChainCount);
	}

void DBSearcher::SetProfiles()
	{
	StartTimer(SetProfiles);
	asserta(m_ThreadCount != UINT_MAX);
	asserta(m_ThreadCount != 0);
	DSS &D = m_D;
	m_Profiles.clear();
	m_ComboLettersVec.clear();
	asserta(m_ChainCount < 100000);
	m_Profiles.resize(m_ChainCount);
	m_ComboLettersVec.resize(m_ChainCount);
	for (uint ChainIndex = 0; ChainIndex < m_ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, m_ChainCount, "Set profiles");
		const PDBChain &Chain = *m_Chains[ChainIndex];
		D.Init(Chain);
		D.GetProfile(m_Profiles[ChainIndex]);
		D.GetComboLetters(m_ComboLettersVec[ChainIndex]);
		}
	EndTimer(SetProfiles);
	}

const PDBChain &DBSearcher::GetChain(uint ChainIndex) const
	{
	asserta(ChainIndex < SIZE(m_Chains));
	return *m_Chains[ChainIndex];
	}

uint DBSearcher::GetQueryCount() const
	{
	return m_QueryChainCount;
	}

bool DBSearcher::GetNextPairQuerySelf(uint &ChainIndex1, uint &ChainIndex2)
	{
	ChainIndex1 = UINT_MAX;
	ChainIndex2 = UINT_MAX;
	m_Lock.lock();
	if (m_PairIndex == m_PairCount)
		{
		m_Lock.unlock();
		return false;
		}
	ProgressStep(m_PairIndex, m_PairCount, "Aligning");
	ChainIndex1 = m_NextChainIndex1;
	ChainIndex2 = m_NextChainIndex2;
	asserta(m_NextChainIndex1 < m_ChainCount);
	asserta(m_NextChainIndex2 < m_ChainCount);
	++m_NextChainIndex2;
	if (m_NextChainIndex2 == m_ChainCount)
		{
		++m_NextChainIndex1;
		m_NextChainIndex2 = m_NextChainIndex1 + 1;
		}
	++m_PairIndex;
	m_Lock.unlock();
	return true;
	}

bool DBSearcher::GetNextPair(uint &QueryChainIndex, uint &DBChainIndex)
	{
	if (m_QuerySelf)
		{
		bool Ok = GetNextPairQuerySelf(QueryChainIndex, DBChainIndex);
		if (Ok)
			asserta(QueryChainIndex != DBChainIndex);
		return Ok;
		}

	QueryChainIndex = UINT_MAX;
	DBChainIndex = UINT_MAX;
	m_Lock.lock();
	if (m_PairIndex == m_PairCount)
		{
		m_Lock.unlock();
		return false;
		}

	asserta(m_PairIndex < m_PairCount);
	ProgressStep(m_PairIndex++, m_PairCount, "Aligning");

	QueryChainIndex = m_NextQueryIdx;
	DBChainIndex = m_QueryChainCount + m_NextDBIdx;
	asserta(QueryChainIndex < m_ChainCount);
	asserta(DBChainIndex < m_ChainCount);

	++m_NextQueryIdx;
	if (m_NextQueryIdx == m_QueryChainCount)
		{
		m_NextQueryIdx = 0;
		++m_NextDBIdx;
		}

	m_Lock.unlock();
	return true;
	}

void DBSearcher::Thread(uint ThreadIndex)
	{
	asserta(ThreadIndex < SIZE(m_DAs));
	uint PrevChainIndex1 = UINT_MAX;
	DSSAligner &DA = *m_DAs[ThreadIndex];
	for (;;)
		{
		uint ChainIndex1, ChainIndex2;
		bool Ok = GetNextPair(ChainIndex1, ChainIndex2);
		if (!Ok)
			break;

		if (ChainIndex1 == PrevChainIndex1)
			++m_QPCacheHits;
		else
			{
			++m_QPCacheMisses;
			const PDBChain &Chain1 = *m_Chains[ChainIndex1];
			const vector<vector<byte> > &Profile1 = m_Profiles[ChainIndex1];
			const vector<byte> &ComboLetters1 = m_ComboLettersVec[ChainIndex1];
			const vector<uint> &KmerBits1 = m_KmerBitsVec[ChainIndex1];
			float Gumbel_mu = FLT_MAX;
			float Gumbel_beta = FLT_MAX;
			if (!m_Gumbel_mus.empty())
				{
				Gumbel_mu = m_Gumbel_mus[ChainIndex1];
				Gumbel_beta = m_Gumbel_betas[ChainIndex1];
				}
			DA.SetQuery(Chain1, &Profile1, &KmerBits1, &ComboLetters1,
			  Gumbel_mu, Gumbel_beta);
			}

		const PDBChain &Chain2 = *m_Chains[ChainIndex2];
		const vector<vector<byte> > &Profile2 = m_Profiles[ChainIndex2];
		const vector<byte> &ComboLetters2 = m_ComboLettersVec[ChainIndex2];
		const vector<uint> &KmerBits2 = m_KmerBitsVec[ChainIndex2];
		float Gumbel_mu = FLT_MAX;
		float Gumbel_beta = FLT_MAX;
		if (!m_Gumbel_mus.empty())
			{
			Gumbel_mu = m_Gumbel_mus[ChainIndex2];
			Gumbel_beta = m_Gumbel_betas[ChainIndex2];
			}
		DA.SetTarget(Chain2, &Profile2, &KmerBits2, &ComboLetters2,
		  Gumbel_mu, Gumbel_beta);

		if (m_Params->m_UseComboPath)
			{
			DA.AlignComboOnly();
			m_Lock.lock();
			OnAln(ChainIndex1, ChainIndex2, DA, true);
			OnAln(ChainIndex1, ChainIndex2, DA, false);
			m_Lock.unlock();
			}
		else if (m_Params->m_ComboScoreOnly)
			{
			float ComboScore = DA.GetComboScore();
			DA.m_EvalueA = ComboScore;
			DA.m_EvalueB = ComboScore;
			m_Lock.lock();
			OnAln(ChainIndex1, ChainIndex2, DA, true);
			OnAln(ChainIndex1, ChainIndex2, DA, false);
			m_Lock.unlock();
			}
		else
			{
			DA.AlignQueryTarget();
			if (!DA.m_Path.empty())
				{
				if (m_CollectTestStats)
					{
					m_Lock.lock();
					asserta(ChainIndex1 < SIZE(m_TestStatsVec));
					vector<float> &v1 = m_TestStatsVec[ChainIndex1];
					v1.push_back(DA.m_TestStatisticA);
					m_Lock.unlock();
					}

				m_Lock.lock();
				DA.ToTsv(m_fTsv, m_MaxEvalue, true);
				DA.ToAln(m_fAln, m_MaxEvalue, true);
				DA.ToFasta2(m_fFasta2, m_MaxEvalue, opt_global, true);
				OnAln(ChainIndex1, ChainIndex2, DA, true);
				m_Lock.unlock();
				if (m_QuerySelf)
					{
					if (m_CollectTestStats)
						{
						m_Lock.lock();
						asserta(ChainIndex2 < SIZE(m_TestStatsVec));
						vector<float> &v2 = m_TestStatsVec[ChainIndex2];
						v2.push_back(DA.m_TestStatisticB);
						m_Lock.unlock();
						}

					m_Lock.lock();
					DA.ToTsv(m_fTsv, m_MaxEvalue, false);
					DA.ToAln(m_fAln, m_MaxEvalue, false);
					OnAln(ChainIndex1, ChainIndex2, DA, false);
					m_Lock.unlock();
					}
				}
			}
		PrevChainIndex1 = ChainIndex1;
		}
	}

uint DBSearcher::GetDBSize() const
	{
	if (m_QuerySelf)
		return m_QueryChainCount;
	else
		return m_DBChainCount;
	}

void DBSearcher::StaticThread(uint ThreadIndex, DBSearcher *ptrDBS)
	{
	ptrDBS->Thread(ThreadIndex);
	}

void DBSearcher::Run()
	{
	m_D.m_Params = m_Params;
	for (uint i = 0; i < SIZE(m_DAs); ++i)
		m_DAs[i]->m_Params = m_Params;

	if (m_Params->m_USort)
		{
		RunUSort();
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
	m_ProcessedQueryCount = UINT_MAX;

	m_PairIndex = 0;
	if (m_QuerySelf)
		{
		asserta(m_ChainCount == m_QueryChainCount);
		asserta(m_DBChainCount == 0);
		m_PairCount = (m_ChainCount*(m_ChainCount-1))/2;
		m_NextChainIndex1 = 0;
		m_NextChainIndex2 = 1;
		}
	else
		{
		m_PairCount = m_QueryChainCount*m_DBChainCount;
		m_NextQueryIdx = 0;
		m_NextDBIdx = 0;
		}

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThread, ThreadIndex, this);
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
	RunStats();
	}

void DBSearcher::RunStats() const
	{
	ProgressLog("Search time %u secs, ", m_Secs);
	DSSAligner::Stats();
	uint Hits = m_QPCacheHits;
	uint Misses = m_QPCacheMisses;
	ProgressLog("QP cache hits %u, misses %u\n", Hits, Misses);
	}

uint DBSearcher::GetDBChainIndex(uint Idx) const
	{
	if (m_QuerySelf)
		return Idx;
	return m_QueryChainCount + Idx;
	}

void DBSearcher::Setup(const DSSParams &Params)
	{
	if (optset_evalue)
		m_MaxEvalue = (float) opt_evalue;
	else
		m_MaxEvalue = FLT_MAX;

	uint ThreadCount = GetRequestedThreadCount();
	asserta(ThreadCount > 0);
	asserta(ThreadCount != UINT_MAX);
	asserta(m_ThreadCount == UINT_MAX);
	asserta(m_DAs.empty());
	asserta(m_Mems.empty());

	m_ProcessedQueryCount = 0;

	m_Params = &Params;
	m_D.m_Params = &Params;
	m_ThreadCount = ThreadCount;

	for (uint i = 0; i < ThreadCount; ++i)
		{
		DSSAligner *DA = new DSSAligner;
		DA->m_Params = m_Params;
		m_DAs.push_back(DA);

		XDPMem *Mem = new XDPMem;
		m_Mems.push_back(Mem);
		}

	vector<FEATURE> ComboFeatures;
	ComboFeatures.push_back(FEATURE_SS3);
	ComboFeatures.push_back(FEATURE_NbrSS3);
	ComboFeatures.push_back(FEATURE_RevNbrDist4);
	DSSParams::SetComboFeatures(ComboFeatures);

	m_D.m_Params = m_Params;
	const uint AS = m_D.GetAlphaSize(FEATURE_Combo);
	m_D.m_PatternAlphaSize1 = AS;
	uint PatternOnes = GetPatternOnes(m_D.m_Params->m_PatternStr);
	m_D.m_PatternAlphaSize = myipow(AS, PatternOnes);
	SetProfiles();
	SetKmersVec();

	if (m_CollectTestStats)
		{
		ProgressLog("\n --- collect teststats ---\n\n");
		m_TestStatsVec.clear();
		m_TestStatsVec.resize(m_ChainCount);
		}
#if SLOPE_CALIB
	LoadCalibratedSlopes(opt_slopes);
#endif
#if GUMBEL_CALIB
	LoadGumbelCalib(opt_gumin);
#endif
	OnSetup();
	}

uint DBSearcher::GetDBChainCount() const
	{
	if (m_QuerySelf)
		return GetQueryChainCount();
	asserta(m_DBChainCount != UINT_MAX);
	asserta(m_DBChainCount > 0);
	return m_DBChainCount;
	}

uint DBSearcher::GetQueryChainCount() const
	{
	asserta(m_QueryChainCount != UINT_MAX);
	asserta(m_QueryChainCount > 0);
	return m_QueryChainCount;
	}

uint DBSearcher::GetQueryChainIdx(uint Idx) const
	{
	asserta(Idx < m_QueryChainCount);
	return Idx;
	}

uint DBSearcher::GetDBChainIdx(uint Idx) const
	{
	if (m_QuerySelf)
		{
		asserta(Idx < m_QueryChainCount);
		return Idx;
		}
	else
		{
		asserta(Idx < m_DBChainCount);
		asserta(m_QueryChainCount + Idx < SIZE(m_Chains));
		return m_QueryChainCount + Idx;
		}
	}

const PDBChain &DBSearcher::GetDBChain(uint Idx) const
	{
	uint ChainIdx = GetDBChainIdx(Idx);
	return GetChain(ChainIdx);
	}

const PDBChain &DBSearcher::GetQueryChain(uint Idx) const
	{
	uint ChainIdx = GetQueryChainIdx(Idx);
	return GetChain(ChainIdx);
	}

const char *DBSearcher::GetQueryLabel(uint Idx) const
	{
	const PDBChain &Chain = GetQueryChain(Idx);
	return Chain.m_Label.c_str();
	}

const char *DBSearcher::GetDBLabel(uint Idx) const
	{
	const PDBChain &Chain = GetDBChain(Idx);
	return Chain.m_Label.c_str();
	}

#if SLOPE_CALIB
void DBSearcher::GetChainSlope(uint ChainIdx, float &m, float &b) const
	{
	if (m_ms.empty())
		{
		m = FLT_MAX;
		b = FLT_MAX;
		return;
		}
	asserta(ChainIdx < SIZE(m_ms));
	asserta(ChainIdx < SIZE(m_bs));
	m = m_ms[ChainIdx];
	b = m_bs[ChainIdx];
	}

void DBSearcher::LoadCalibratedSlopes(const string &FN)
	{
	m_ms.clear();
	m_bs.clear();
	if (FN == "")
		return;
	m_ms.resize(m_ChainCount, FLT_MAX);
	m_bs.resize(m_ChainCount, FLT_MAX);
	FILE *f = OpenStdioFile(FN);
	string Line;
	vector<string> Fields;
	map<string, uint> LabelToIdx;
	vector<float> ms;
	vector<float> bs;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		if (Fields[1] == "m")
			continue;
		asserta(SIZE(Fields) == 3);
		const string &Label = Fields[0];
		LabelToIdx[Label] = SIZE(ms);
		float m = (float) StrToFloat(Fields[1]);
		float b = (float) StrToFloat(Fields[2]);
		ms.push_back(m);
		bs.push_back(b);
		}
	CloseStdioFile(f);
	for (uint ChainIdx = 0; ChainIdx < m_ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *m_Chains[ChainIdx];
		const string &Label = Chain.m_Label;
		string Dom;
		SCOP40Bench::GetDomFromLabel(Label, Dom);
		map<string, uint>::iterator iter = LabelToIdx.find(Dom);
		if (iter == LabelToIdx.end())
			Die("Label not found in slopes >%s", Label.c_str());
		uint Idx = iter->second;
		asserta(Idx < SIZE(ms));
		asserta(Idx < SIZE(bs));
		float m = ms[Idx];
		float b = bs[Idx];
		m_ms[ChainIdx] = m;
		m_bs[ChainIdx] = b;
		}
	}
#endif // SLOPE_CALIB

#if GUMBEL_CALIB
void DBSearcher::LoadGumbelCalib(const string &FN)
	{
	m_Gumbel_mus.clear();
	m_Gumbel_betas.clear();
	if (FN == "")
		return;
	m_Gumbel_mus.resize(m_ChainCount, FLT_MAX);
	m_Gumbel_betas.resize(m_ChainCount, FLT_MAX);
	FILE *f = OpenStdioFile(FN);
	string Line;
	vector<string> Fields;
	map<string, uint> LabelToIdx;
	vector<float> mus;
	vector<float> betas;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		if (Fields[1] == "m")
			continue;
		asserta(SIZE(Fields) == 3);
		const string &Label = Fields[0];
		LabelToIdx[Label] = SIZE(mus);
		float mu = (float) StrToFloat(Fields[1]);
		float beta = (float) StrToFloat(Fields[2]);
		mus.push_back(mu);
		betas.push_back(beta);
		}
	CloseStdioFile(f);
	for (uint ChainIdx = 0; ChainIdx < m_ChainCount; ++ChainIdx)
		{
		const PDBChain &Chain = *m_Chains[ChainIdx];
		const string &Label = Chain.m_Label;
		map<string, uint>::iterator iter = LabelToIdx.find(Label);
		if (iter == LabelToIdx.end())
			Die("Label not found in gumin >%s", Label.c_str());
		uint Idx = iter->second;
		asserta(Idx < SIZE(mus));
		asserta(Idx < SIZE(betas));
		float mu = mus[Idx];
		float beta = betas[Idx];
		m_Gumbel_mus[ChainIdx] = mu;
		m_Gumbel_betas[ChainIdx] = beta;
		}
	}

void DBSearcher::GetChainGumbel(uint ChainIdx, float &mu, float &beta) const
	{
	if (m_Gumbel_mus.empty())
		{
		mu = FLT_MAX;
		beta = FLT_MAX;
		return;
		}
	asserta(ChainIdx < SIZE(m_Gumbel_mus));
	asserta(ChainIdx < SIZE(m_Gumbel_betas));
	mu = m_Gumbel_mus[ChainIdx];
	beta = m_Gumbel_betas[ChainIdx];
	}
#endif
