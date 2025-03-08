#include "myutils.h"
#include "dssparams.h"
#include "muprefilter.h"
#include "muprefilter_params.h"
#include "museqsource.h"
#include "seqinfo.h"
#include "mymutex.h"

#define USE_SEQ_FEEDER 0

#if USE_SEQ_FEEDER
#include "seqfeeder.h"
static SeqFeeder *s_SF;
#endif

static const MerMx *s_ptrScoreMx;
static const SeqDB *s_ptrQDB = 0;
static const MuDex *s_ptrQKmerIndex = 0;
static time_t s_TimeLastProgress;
static uint s_DBSize = 0;
static MuSeqSource *s_SS;
static mutex s_ProgressLock;

static void ThreadBody_Filter(uint ThreadIndex)
	{
	MuPrefilter Pref;
	Pref.m_ScoreMx = s_ptrScoreMx;
	Pref.m_QKmerIndex = s_ptrQKmerIndex;
	Pref.m_KmerSelfScores = s_ptrQKmerIndex->m_KmerSelfScores;
	Pref.SetQDB(*s_ptrQDB);
#if !USE_SEQ_FEEDER
	ObjMgr OM;
#endif

	for (;;)
		{
#if USE_SEQ_FEEDER
		SeqInfo *SI = s_SF->GetSI(ThreadIndex);
		if (SI == 0)
			return;
#else
		SeqInfo *SI = OM.GetSeqInfo();
		bool Ok = s_SS->GetNext(SI);
		if (!Ok)
			return;
#endif

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
#if USE_SEQ_FEEDER
			s_SF->Down(ThreadIndex, SI);
#else
			OM.Down(SI);
#endif
			continue;
			}

		const byte *TSeq = SI->m_Seq;
		const string &TLabel = string(SI->m_Label);
		uint TL = SI->m_L;
		Pref.Search(0, TSeqIdx, TLabel, TSeq, TL);
#if USE_SEQ_FEEDER
		s_SF->Down(ThreadIndex, SI);
#else
		OM.Down(SI);
#endif
		}
	}

uint MuPreFilter(const DSSParams &Params,
			  SeqDB &QDB,
			  MuSeqSource &FSS,
			  const string &OutputFN)
	{
	s_SS = &FSS;
	s_SS->m_ASCII = false;

	const uint k = MuDex::m_k;

	QDB.ToLetters(g_CharToLetterMu);
	const uint QSeqCount = QDB.GetSeqCount();

	MuPrefilter::m_RSB.m_B = RSB_SIZE;
	MuPrefilter::m_RSB.Init(QSeqCount);

	const MerMx &ScoreMx = GetMuMerMx(k);
	asserta(ScoreMx.m_k == k);

	MuDex QKmerIndex;
	QKmerIndex.m_KmerSelfScores = ScoreMx.BuildSelfScores_Kmers();
	QKmerIndex.m_MinKmerSelfScore = MIN_KMER_PAIR_SCORE;
	QKmerIndex.FromSeqDB(QDB);
#if DEBUG
	QKmerIndex.Validate();
#endif
	asserta(QKmerIndex.m_k == k);
	asserta(QKmerIndex.m_DictSize == DICT_SIZE);
	asserta(ScoreMx.m_AS_pow[k] == QKmerIndex.m_DictSize);

	s_ptrQDB = &QDB;
	s_ptrScoreMx = &ScoreMx;
	s_ptrQKmerIndex = &QKmerIndex;

	Progress("Filtering\r");
	time_t t_start = time(0);
	s_TimeLastProgress = t_start;

	uint ThreadCount = GetRequestedThreadCount();
#if USE_SEQ_FEEDER
	s_SF = new SeqFeeder;
	s_SF->m_SS = &FSS;
	s_SF->Start(ThreadCount);
#else
	s_SS = &FSS;
#endif

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
	Progress("Filtering done %s      \n", IntToStr(s_DBSize));
#if SEQ_FEEDER_STATS
	s_SF->Stats();
#endif

	FILE *fTsv = CreateStdioFile(OutputFN);
	MuPrefilter::m_RSB.ToTsv(fTsv);
	CloseStdioFile(fTsv);
	return s_DBSize;
	}
