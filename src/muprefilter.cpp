#include "myutils.h"
#include "dssparams.h"
#include "prefiltermu.h"
#include "prefiltermuparams.h"
#include "museqsource.h"
#include "seqinfo.h"

static mutex s_NextLock;
static const MerMx *s_ptrScoreMx;
static const SeqDB *s_ptrQDB = 0;
static const MuDex *s_ptrQKmerIndex = 0;
static time_t s_TimeLastProgress;
static MuSeqSource *s_SS;
static uint s_DBSize = 0;

#if USE_BIAS
static const vector<vector<int8_t> > *s_ptrBiasVecs8;
#endif

static void ThreadBody(uint ThreadIndex)
	{
	ObjMgr OM;
	PrefilterMu Pref;
	Pref.m_ScoreMx = s_ptrScoreMx;
	Pref.m_QKmerIndex = s_ptrQKmerIndex;
	Pref.m_KmerSelfScores = s_ptrQKmerIndex->m_KmerSelfScores;
	Pref.SetQDB(*s_ptrQDB);
#if USE_BIAS
	Pref.m_BiasVecs8 = s_ptrBiasVecs8;
#endif

	for (;;)
		{
		s_NextLock.lock();
		uint TSeqIdx = s_DBSize;
		SeqInfo *SI = OM.GetSeqInfo();
		bool Ok = s_SS->GetNext(SI);
		if (Ok)
			++s_DBSize;
		if (s_DBSize%100 == 0)
			{
			time_t now = time(0);
			if (now > s_TimeLastProgress)
				Progress("Filtering %s    \r", IntToStr(s_DBSize));
			s_TimeLastProgress = now;
			}
		s_NextLock.unlock();
		if (!Ok)
			{
			OM.Down(SI);
			return;
			}
		if (SI->m_L == 0)
			{
			OM.Down(SI);
			continue;
			}

		const byte *TSeq = SI->m_Seq;
		const string &TLabel = string(SI->m_Label);
		uint TL = SI->m_L;
		Pref.Search(0, TSeqIdx, TLabel, TSeq, TL);
		OM.Down(SI);
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

	PrefilterMu::m_RSB.m_B = RSB_SIZE;
	PrefilterMu::m_RSB.Init(QSeqCount);

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

#if USE_BIAS
	vector<vector<int8_t> > BiasVecs8(QSeqCount);
	vector<float> BiasVec;
	for (uint QSeqIdx = 0; QSeqIdx < QSeqCount; ++QSeqIdx)
		{
		vector<int8_t> &BiasVec8 = BiasVecs8[QSeqIdx];
		const byte *QSeq = QDB.GetByteSeq(QSeqIdx);
		uint QL = QDB.GetSeqLength(QSeqIdx);
		CalcLocalBiasCorrection_Mu(QSeq, QL, BIAS_WINDOW, QBIAS_SCALE, BiasVec, BiasVec8);
		}
#endif

	s_ptrQDB = &QDB;
	s_ptrScoreMx = &ScoreMx;
	s_ptrQKmerIndex = &QKmerIndex;
#if USE_BIAS
	s_ptrBiasVecs8 = &BiasVecs8;
#endif

	Progress("Filtering\r");
	time_t t_start = time(0);
	s_TimeLastProgress = t_start;

	vector<thread *> ts;
	uint ThreadCount = GetRequestedThreadCount();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(ThreadBody, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	Progress("Filtering done       \n");

	FILE *fTsv = CreateStdioFile(OutputFN);
	PrefilterMu::m_RSB.ToTsv(fTsv);
	CloseStdioFile(fTsv);
	return s_DBSize;
	}
