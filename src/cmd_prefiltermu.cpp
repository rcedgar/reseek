#include "myutils.h"
#include "prefiltermu.h"
#include "prefiltermuparams.h"

static uint s_NextTIdx = 0;
static mutex m_NextTIdxLock;
static const MerMx *s_ptrScoreMx;
static const SeqDB *s_ptrQDB = 0;
static const SeqDB *s_ptrTDB = 0;
static const vector<vector<int8_t> > *s_ptrBiasVecs8;
static const MuDex *s_ptrQKmerIndex = 0;
static FILE *s_fTsv = 0;
static FILE *s_fTsv2 = 0;
static time_t s_TimeLastProgress;

static void ThreadBody(uint ThreadIndex)
	{
	const uint TSeqCount = s_ptrTDB->GetSeqCount();

	PrefilterMu Pref;
	Pref.m_ScoreMx = s_ptrScoreMx;
	Pref.m_QKmerIndex = s_ptrQKmerIndex;
	Pref.m_KmerSelfScores = s_ptrQKmerIndex->m_KmerSelfScores;
	Pref.SetQDB(*s_ptrQDB);
	Pref.m_BiasVecs8 = s_ptrBiasVecs8;

	for (;;)
		{
		m_NextTIdxLock.lock();
		uint TSeqIdx = s_NextTIdx;
		if (s_NextTIdx < TSeqCount)
			++s_NextTIdx;
		if (TSeqIdx > 0 && TSeqIdx + 1 < TSeqCount)
			{
			time_t now = time(0);
			if (now > s_TimeLastProgress)
				ProgressStep(TSeqIdx, TSeqCount, "Filtering");
			s_TimeLastProgress = now;
			}
		m_NextTIdxLock.unlock();
		if (TSeqIdx == TSeqCount)
			return;

		Pref.m_TSeqIdx = TSeqIdx;
		const byte *TSeq = s_ptrTDB->GetByteSeq(TSeqIdx);
		const string &TLabel = s_ptrTDB->GetLabel(TSeqIdx);
		uint TL = s_ptrTDB->GetSeqLength(TSeqIdx);
		Pref.Search(s_fTsv2, TSeqIdx, TLabel, TSeq, TL);
		}
	}

void cmd_prefilter_mu()
	{
	const uint k = MuDex::m_k;

	const string &Query3Di_FN = g_Arg1;
	const string &DB3Di_FN = opt_db;

	s_fTsv2 = CreateStdioFile(opt_output2);
	
	SeqDB QDB;
	SeqDB TDB;

	QDB.FromFasta(Query3Di_FN);
	TDB.FromFasta(DB3Di_FN);

	QDB.ToLetters(g_CharToLetterMu);
	TDB.ToLetters(g_CharToLetterMu);
	const uint QSeqCount = QDB.GetSeqCount();
	const uint TSeqCount = TDB.GetSeqCount();

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

	vector<vector<int8_t> > BiasVecs8(QSeqCount);
	vector<float> BiasVec;
	for (uint QSeqIdx = 0; QSeqIdx < QSeqCount; ++QSeqIdx)
		{
		vector<int8_t> &BiasVec8 = BiasVecs8[QSeqIdx];
		const byte *QSeq = QDB.GetByteSeq(QSeqIdx);
		uint QL = QDB.GetSeqLength(QSeqIdx);
		CalcLocalBiasCorrection_Mu(QSeq, QL, BIAS_WINDOW, QBIAS_SCALE, BiasVec, BiasVec8);
		}

	s_ptrQDB = &QDB;
	s_ptrTDB = &TDB;
	s_ptrScoreMx = &ScoreMx;
	s_ptrQKmerIndex = &QKmerIndex;
	s_ptrBiasVecs8 = &BiasVecs8;

	ProgressStep(0, TSeqCount, "Filtering");
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
	ProgressStep(TSeqCount-1, TSeqCount, "Filtering");
	CloseStdioFile(s_fTsv2);

	time_t t_end = time(0);
	uint filter_secs = uint(t_end - t_start);
	double SeqsPerSec = double(TSeqCount)/filter_secs;
	double SeqsPerSecPerThread = SeqsPerSec*ThreadCount;
	ProgressLog("Seqs/sec         %s\n", FloatToStr(SeqsPerSec));
	ProgressLog("Seqs/sec/thread  %s\n", FloatToStr(SeqsPerSecPerThread));

	{
	FILE *fTsv = CreateStdioFile(opt_output);
	PrefilterMu::m_RSB.ToTsv(fTsv);
	CloseStdioFile(s_fTsv);
	}

	if (optset_output3)
		{
		vector<string> QLabels;
		vector<string> TLabels;
		for (uint i = 0; i < QSeqCount; ++i)
			QLabels.push_back(QDB.GetLabel(i));
		for (uint i = 0; i < TSeqCount; ++i)
			TLabels.push_back(TDB.GetLabel(i));
		FILE *fTsv = CreateStdioFile(opt_output3);
		PrefilterMu::m_RSB.ToLabelsTsv(fTsv, QLabels, TLabels);
		CloseStdioFile(s_fTsv);
		}
	}
