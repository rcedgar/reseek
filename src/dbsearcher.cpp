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
			m_MaxEvalue = 1e6;
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

	OnSetup();
	}

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

void DBSearcher::BaseOnAln(DSSAligner &DA, bool Up)
	{
	if (!opt(scores_are_not_evalues) && DA.GetEvalue(Up) > m_MaxEvalue)
		return;
	m_Lock.lock();
	++m_HitCount;
	DA.ToTsv(g_fTsv, Up);
	DA.ToAln(g_fAln, Up);
	DA.ToFasta2(g_fFasta2, Up);
	OnAln(DA, Up);
	m_Lock.unlock();
	}
