#include "myutils.h"
#include "cigar.h"
#include "dss.h"
#include "dbsearcher.h"
#include "profileloader.h"
#include "chainreader2.h"
#include "scop40bench.h"
#include <thread>
#include "timing.h"

#if 0
void DBSearcher::SetKmersVec()
	{
	StartTimer(SetKmersVec);
	DSS &D = m_D;
	m_KmersVec.clear();
	m_KmerBitsVec.clear();
	asserta(m_ChainCount < 100000);
	m_KmersVec.resize(m_ChainCount);
	m_KmerBitsVec.resize(m_ChainCount);
	vector<byte> ComboLetters;
	for (uint ChainIndex = 0; ChainIndex < m_ChainCount; ++ChainIndex)
		{
		ProgressStep(ChainIndex, m_ChainCount, "Set k-mers");
		const PDBChain &Chain = *m_Chains[ChainIndex];
		D.Init(Chain);
		vector<uint> *Kmers = new vector<uint>;
		vector<uint> *Bits = new vector<uint>;
		D.GetComboLetters(ComboLetters);
		D.GetComboKmers(ComboLetters, *Kmers);
		D.GetComboKmerBits(*Kmers, *Bits);
		m_KmersVec[ChainIndex] = Kmers;
		m_KmerBitsVec[ChainIndex] = Bits;
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
		vector<vector<byte> > *ptrProfile = new vector<vector<byte> >;
		vector<byte> *ptrComboLetters = new vector<byte>;
		D.GetProfile(*ptrProfile);
		m_Profiles[ChainIndex] = ptrProfile;
		D.GetComboLetters(*ptrComboLetters);
		m_ComboLettersVec[ChainIndex] = ptrComboLetters;
		}
	EndTimer(SetProfiles);
	}
#endif // 0

uint DBSearcher::GetDBSize() const
	{
	return GetDBChainCount();
	}


void DBSearcher::RunStats() const
	{
	uint Secs = m_Secs;
	if (Secs == 0)
		Secs = 1;
	uint ThreadCount = GetRequestedThreadCount();
	double QueriesPerSec = m_ProcessedQueryCount/Secs;
	double PairsPerSec = m_ProcessedPairCount/Secs;
	double PairsPerSecPerThread = PairsPerSec/ThreadCount;
	ProgressLog("\n");
	ProgressLog("%10.10s  Search time\n", SecsToHHMMSS(Secs));
	ProgressLog("%10.10s  Query chains\n", IntToStr(m_ProcessedQueryCount));
	ProgressLog("%10.10s  DB chains\n", IntToStr(GetDBChainCount()));
	ProgressLog("%10.1f  Queries/sec\n", QueriesPerSec);
	ProgressLog("%10.10s  Query-target comparisons/sec\n", FloatToStr(PairsPerSec));
	ProgressLog("%10.10s  Query-target comparisons/sec/thread (%u threads)\n", FloatToStr(PairsPerSecPerThread), ThreadCount);
	ProgressLog("%10.10s  Alignments\n", IntToStr(DSSAligner::m_AlnCount));
	ProgressLog("%10.10s  Hits (max E-value %.3g)\n", IntToStr(m_HitCount), m_MaxEvalue);

	DSSAligner::Stats();
	uint Hits = m_QPCacheHits;
	uint Misses = m_QPCacheMisses;
	Log("QP cache hits %u, misses %u\n", Hits, Misses);
	}

void DBSearcher::Setup()
	{
	if (optset_evalue)
		m_MaxEvalue = (float) opt_evalue;
	else
		m_MaxEvalue = 10;

	uint ThreadCount = GetRequestedThreadCount();
	asserta(ThreadCount > 0);
	asserta(ThreadCount != UINT_MAX);
	asserta(m_ThreadCount == UINT_MAX);
	asserta(m_DAs.empty());
	asserta(m_Mems.empty());

	m_ProcessedQueryCount = 0;
	m_ProcessedPairCount = 0;
	m_HitCount = 0;

	m_D.SetParams(*m_Params);
	m_ThreadCount = ThreadCount;

	for (uint i = 0; i < ThreadCount; ++i)
		{
		DSSAligner *DA = new DSSAligner;
		DA->m_Params = m_Params;
		m_DAs.push_back(DA);

		XDPMem *Mem = new XDPMem;
		m_Mems.push_back(Mem);
		}
#if SLOPE_CALIB
	LoadCalibratedSlopes(opt_slopes);
#endif
#if GUMBEL_CALIB
	LoadGumbelCalib(opt_gumin);
#endif
	OnSetup();
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
	uint ChainCount = GetDBChainCount();
	m_Gumbel_mus.resize(ChainCount, FLT_MAX);
	m_Gumbel_betas.resize(ChainCount, FLT_MAX);
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
		const PDBChain &Chain = *m_DBChains[ChainIdx];
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

void DBSearcher::LoadDB(const string &DBFN)
	{
	ChainReader2 CR;
	CR.Open(DBFN);
	ProfileLoader PL;
	uint ThreadCount = GetRequestedThreadCount();

	vector<PDBChain *> *ptrChains = &m_DBChains;
	vector<vector<vector<byte> > *> *ptrProfiles = &m_DBProfiles;
	vector<vector<byte> *> *ptrComboLetters = &m_DBComboLettersVec;
	vector<vector<uint> *> *ptrKmerBitsVec = &m_DBKmerBitsVec;

	if (m_Params->m_MinU <= 0)
		ptrKmerBitsVec = 0;
	if (m_Params->m_Omega <= 0)
		ptrComboLetters = 0;

	PL.Load(CR, ptrChains, ptrProfiles, ptrComboLetters,
	  ptrKmerBitsVec, *m_Params, ThreadCount);
	}
