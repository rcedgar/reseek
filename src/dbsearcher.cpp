#include "myutils.h"
#include "cigar.h"
#include "dss.h"
#include "dbsearcher.h"
#include "profileloader.h"
#include "chainreader2.h"
#include "scop40bench.h"
#include <thread>
#include "timing.h"
#include "output.h"

DBSearcher::~DBSearcher()
	{
#define delv(x)	for (uint i = 0; i < SIZE(x); ++i) delete x[i];

	delv(m_Mems)
	delv(m_DAs)
	delv(m_DBChains)
	delv(m_DBProfiles);
	delv(m_DBMuLettersVec);
	delv(m_DBMuKmersVec);
	}

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
	if (m_MaxEvalue == DBL_MAX)
		ProgressLog("%10.10s  Hits\n", IntToStr(m_HitCount), m_MaxEvalue);
	else
		ProgressLog("%10.10s  Hits (max E-value %.3g)\n", IntToStr(m_HitCount), m_MaxEvalue);
	if (m_ProcessedQueryCount < 100)
		return;
	ProgressLog("%10.10s  DB chains\n", IntToStr(m_ProcessedQueryCount));
	ProgressLog("%10.10s  Query chains\n", IntToStr(GetDBChainCount()));
	ProgressLog("%10.1f  Chains/sec\n", QueriesPerSec);
	ProgressLog("%10.10s  Comparisons/sec\n", FloatToStr(PairsPerSec));
	ProgressLog("%10.10s  Comparisons/sec/thread (%u threads)\n", FloatToStr(PairsPerSecPerThread), ThreadCount);

	DSSAligner::Stats();
	uint Hits = m_QPCacheHits;
	uint Misses = m_QPCacheMisses;
	Log("QP cache hits %u, misses %u\n", Hits, Misses);
	}

void DBSearcher::AddChain(PDBChain *ptrChain, vector<vector<byte> > *ptrProfile,
  vector<byte> *ptrMuLetters)
	{
	m_DBChains.push_back(ptrChain);
	m_DBProfiles.push_back(ptrProfile);
	m_DBMuLettersVec.push_back(ptrMuLetters);
	}

void DBSearcher::InitEmpty()
	{
// Initialized by c'tor, nothing to do here because
//   we assume object is never re-used
	asserta(SIZE(m_DBChains) == 0);
	}

void DBSearcher::Setup()
	{
	if (optset_evalue)
		m_MaxEvalue = (float) opt(evalue);
	else 
		{
		if (optset_verysensitive)
			m_MaxEvalue = DBL_MAX;
		else
			m_MaxEvalue = 10;
		}

	uint ThreadCount = GetRequestedThreadCount();
	asserta(ThreadCount > 0);
	asserta(ThreadCount != UINT_MAX);
	asserta(m_ThreadCount == UINT_MAX);
	asserta(m_DAs.empty());
	asserta(m_Mems.empty());

	m_ProcessedQueryCount = 0;
	m_ProcessedPairCount = 0;
	m_HitCount = 0;

	m_ThreadCount = ThreadCount;

	for (uint i = 0; i < ThreadCount; ++i)
		{
		DSSAligner *DA = new DSSAligner;
		DA->SetParams(*m_Params);
		m_DAs.push_back(DA);

		XDPMem *Mem = new XDPMem;
		m_Mems.push_back(Mem);
		}

#if SLOPE_CALIB
	LoadCalibratedSlopes(opt(slopes));
#endif
#if GUMBEL_CALIB
	LoadGumbelCalib(opt(gumin));
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

	vector<vector<byte> *> *ptrMuLetters = &m_DBMuLettersVec;
	vector<vector<uint> *> *ptrMuKmersVec = &m_DBMuKmersVec;
	if (m_Params->m_Omega <= 0)
		ptrMuLetters = 0;

	PL.Load(*m_Params, CR, &m_DBChains, &m_DBProfiles, ptrMuLetters,
	  ptrMuKmersVec, &m_DBSelfRevScores, ThreadCount);
	}

bool DBSearcher::Reject(DSSAligner &DA, bool Up) const
	{
	bool Evalue_ok = true;
	bool TS_ok = true;
	if (!opt(scores_are_not_evalues) && DA.GetEvalue(Up) > m_MaxEvalue)
		Evalue_ok = false;
	if (optset_mints && DA.GetNewTestStatistic(Up) < opt(mints))
		TS_ok = false;
	if (Evalue_ok || TS_ok)
		return false;
	return true;
	}

void DBSearcher::BaseOnAln(DSSAligner &DA, bool Up)
	{
	if (Reject(DA, Up))
		return;
	m_Lock.lock();
	++m_HitCount;
	DA.ToTsv(g_fTsv, Up);
	DA.ToAln(g_fAln, Up);
	DA.ToFasta2(g_fFasta2, opt(unaligned), Up);
	OnAln(DA, Up);
	m_Lock.unlock();
	}
