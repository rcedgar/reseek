#include "myutils.h"
#include "dssparams.h"
#include "prefiltermu.h"
#include "prefiltermuparams.h"
#include "museqsource.h"
#include "seqinfo.h"
#include "mymutex.h"

void decide_query_or_db_kmer_neighborhood(uint QSeqCount, uint TSeqCount);

static const MerMx *s_ptrScoreMx;
static const SeqDB *s_ptrQDB = 0;
static const MuDex *s_ptrQKmerIndex = 0;
static time_t s_TimeLastProgress;
static uint s_DBSize = 0;
static MuSeqSource *s_SS;
static mutex s_ProgressLock;

// If true, index k-mers in query+neighborhood
// If false, construct target k-mer neighborhoods.
bool g_QueryNeighborhood = true;

static void ThreadBody_Filter(uint ThreadIndex)
	{
	PrefilterMu Pref;
	Pref.m_OneHitDiag = opt(onehitdiag);
	Pref.m_ScoreMx = s_ptrScoreMx;
	Pref.m_QKmerIndex = s_ptrQKmerIndex;
	Pref.m_KmerSelfScores = s_ptrQKmerIndex->m_KmerSelfScores;
	Pref.SetQDB(*s_ptrQDB);
	ObjMgr OM;

	for (;;)
		{
		SeqInfo *SI = OM.GetSeqInfo();
		bool Ok = s_SS->GetNext(SI);
		if (!Ok)
			return;

		uint TSeqIdx = s_DBSize;
		s_ProgressLock.lock();
		if (s_DBSize%100 == 0)
			{
			time_t now = time(0);
			if (now > s_TimeLastProgress)
				Progress("Filtering %s    \r", IntToStr(s_DBSize));
			s_TimeLastProgress = now;
			}
		++s_DBSize;
		s_ProgressLock.unlock();

		if (SI->m_L == 0)
			{
			OM.Down(SI);
			continue;
			}

		const byte *TSeq = SI->m_Seq;
		const string &TLabel = string(SI->m_Label);
		uint TL = SI->m_L;
		Pref.Search(TSeqIdx, TLabel, TSeq, TL);
		OM.Down(SI);
		}
	}

void decide_query_or_db_kmer_neighborhood(uint QSeqCount, uint DBSeqCount)
	{
	if (opt(idxq))
		g_QueryNeighborhood = true;
	else if (opt(idxt))
		g_QueryNeighborhood = false;
	else
		{
		if (QSeqCount <= MAX_QUERY_CHAINS_FOR_QUERY_NEIGHBORHOOD)
			g_QueryNeighborhood = true;
		else
			g_QueryNeighborhood = false;
		}
	Log("g_QueryNeighborhood=%c\n", tof(g_QueryNeighborhood));
	}

void MuPreFilter(SeqDB &QDB, MuSeqSource &FSS, vector<uint> &TargetIdxs,
	map<uint, vector<uint> > &TargetIdxToQueryIdxs)
	{
	const uint QSeqCount = QDB.GetSeqCount();
	decide_query_or_db_kmer_neighborhood(QSeqCount, UINT_MAX);

	s_SS = &FSS;
	s_SS->m_ASCII = false;

	const uint k = MuDex::m_k;

	QDB.ToLetters(g_CharToLetterMu);
#if DEBUG_CHECKS
	{
	for (uint SeqIdx = 0; SeqIdx < QSeqCount; ++SeqIdx)
		{
		const uint L = QDB.GetSeqLength(SeqIdx);
		const byte *MuLetters = QDB.GetByteSeq(SeqIdx);
		for (uint Pos = 0; Pos < L; ++Pos)
			asserta(MuLetters[Pos] < 36);
		}
	}
#endif

	PrefilterMu::m_RSB.m_B = DSSParams::m_rsb_size;
	PrefilterMu::m_RSB.Init(QSeqCount);

	const MerMx &ScoreMx = GetMuMerMx(k);
	asserta(ScoreMx.m_k == k);

	MuDex QKmerIndex;
	QKmerIndex.m_AddNeighborhood = g_QueryNeighborhood;
	QKmerIndex.m_KmerSelfScores = ScoreMx.BuildSelfScores_Kmers();
	QKmerIndex.m_MinKmerSelfScore =  DSSParams::m_PrefilterMinMuKmerPairScore;
	QKmerIndex.FromSeqDB(QDB);
	asserta(QKmerIndex.m_k == k);
	asserta(QKmerIndex.m_DictSize == PREFILTER_KMER_DICT_SIZE);
	asserta(ScoreMx.m_AS_pow[k] == QKmerIndex.m_DictSize);

	s_ptrQDB = &QDB;
	s_ptrScoreMx = &ScoreMx;
	s_ptrQKmerIndex = &QKmerIndex;

	Progress("Filtering\r");
	time_t t_start = time(0);
	s_TimeLastProgress = t_start;

	uint ThreadCount = GetRequestedThreadCount();
	s_SS = &FSS;

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody_Filter, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	uint Total = PrefilterMu::m_RSB.TruncateAllQueryVecs();
	if (optset_output4)
		{
		const string fn = opt(output4);
		FILE *f = CreateStdioFile(fn);
		ProgressLog("Writing prefilter hits to %s\n", fn.c_str());
		PrefilterMu::m_RSB.ToTsv(f);
		CloseStdioFile(f);
		}
	ProgressLog("%s prefilter hits\n", FloatToStr(double(Total)));
	PrefilterMu::m_RSB.GetTargetInfo(TargetIdxs, TargetIdxToQueryIdxs);
	}
