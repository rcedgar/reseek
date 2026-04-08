#include "myutils.h"
#include "dssparams.h"
#include "prefilter_kappa.h"
#include "kappa_mermx.h"
#include "kappa_dex.h"
#include "kappa_prefilter_params.h"
#include <chrono>

static uint s_NextTIdx = 0;
static mutex m_NextTIdxLock;
static const kappa_mermx *s_ptrScoreMx;
static const SeqDB *s_ptrQDB = 0;
static const SeqDB *s_ptrTDB = 0;
static const kappa_dex *s_ptrQKmerIndex = 0;
static FILE *s_fTsv = 0;
static time_t s_TimeLastProgress;

static uint get_pattern_ones(const string &Str)
	{
	uint n = 0;
	for (uint i = 0; i < SIZE(Str); ++i)
		{
		char c = Str[i];
		asserta(c == '0' || c == '1');
		if (c == '1')
			++n;
		}
	return n;
	}

static void ThreadBody(uint ThreadIndex)
	{
	const uint TSeqCount = s_ptrTDB->GetSeqCount();

	prefilter_kappa Pref;
	Pref.m_ScoreMx = s_ptrScoreMx;
	Pref.m_QKmerIndex = s_ptrQKmerIndex;
	Pref.m_KmerSelfScores = s_ptrQKmerIndex->m_KmerSelfScores;
	Pref.SetQDB(*s_ptrQDB);

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
		Pref.Search(TSeqIdx, TLabel, TSeq, TL);
		}
	}

void cmd_prefilter_kappa()
	{
	const string &QueryKappa_FN = g_Arg1;
	const string &DB3Di_FN = opt(db);

	SeqDB QDB;
	SeqDB TDB;

	QDB.FromFasta(QueryKappa_FN);
	TDB.FromFasta(DB3Di_FN);

	QDB.ToLetters(g_CharToLetterMu);
	TDB.ToLetters(g_CharToLetterMu);
	const uint QSeqCount = QDB.GetSeqCount();
	const uint TSeqCount = TDB.GetSeqCount();

	void SetQueryNeighborhood(uint QSeqCount);
	SetQueryNeighborhood(QSeqCount);

	prefilter_kappa::m_RSB.m_B = DSSParams::m_rsb_size;
	prefilter_kappa::m_RSB.Init(QSeqCount);

	if (optset_kappa_kmer_pattern)
		{
		const string s = opt(kappa_kmer_pattern);
		uint k = get_pattern_ones(s);
		uint K = uint(s.size());

		DSSParams::m_PrefilterKappaPattern = s;
		DSSParams::m_PrefilterKappaKmerNrOnes = k; 
		DSSParams::m_PrefilterKappaKmerWidth = K;
		DSSParams::m_PrefilterKappaDictSize = myipow(32, k);
		}

	kappa_dex QKmerIndex;
	QKmerIndex.Init();
	const uint k = kappa_dex::m_k;

	const kappa_mermx &GetKappaMerMx(uint k);
	const kappa_mermx &ScoreMx = GetKappaMerMx(k);
	asserta(ScoreMx.m_k == k);

	QKmerIndex.m_KmerSelfScores = ScoreMx.BuildSelfScores_Kmers();
	QKmerIndex.m_MinKmerSelfScore =  DSSParams::m_PrefilterMinKappaKmerPairScore;
	QKmerIndex.FromSeqDB(QDB);
#if DEBUG
	QKmerIndex.Validate();
#endif
	asserta(QKmerIndex.m_k == k);
	asserta(QKmerIndex.m_DictSize == DSSParams::m_PrefilterKappaDictSize);
	asserta(ScoreMx.m_AS_pow[k] == QKmerIndex.m_DictSize);

	s_ptrQDB = &QDB;
	s_ptrTDB = &TDB;
	s_ptrScoreMx = &ScoreMx;
	s_ptrQKmerIndex = &QKmerIndex;

	ProgressStep(0, TSeqCount, "Filtering");
	time_t t_start = time(0);
	s_TimeLastProgress = t_start;
	auto chrono_start = std::chrono::high_resolution_clock::now();

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

	time_t t_end = time(0);
	uint filter_secs = uint(t_end - t_start);
	auto chrono_end = std::chrono::high_resolution_clock::now();

	double elapsed_ms = std::chrono::duration<double, std::milli>
		(chrono_end - chrono_start).count();
	double SeqsPerMs= double(TSeqCount)/elapsed_ms;
	ProgressLog("Seqs/ms         %s\n", FloatToStr(SeqsPerMs));

	{
	FILE *fTsv = CreateStdioFile(opt(output));
	prefilter_kappa::m_RSB.ToTsv(fTsv);
	CloseStdioFile(s_fTsv);
	}

#if STORE_PAIR_SCORES
	if (optset_output2)
		{
		vector<string> QLabels;
		vector<string> TLabels;
		for (uint i = 0; i < QSeqCount; ++i)
			QLabels.push_back(QDB.GetLabel(i));
		for (uint i = 0; i < TSeqCount; ++i)
			TLabels.push_back(TDB.GetLabel(i));

		FILE *f = CreateStdioFile(opt(output2));
		const vector<vector<uint16_t> > &QueryIdxToTopScoreVec =
			prefilter_kappa::m_RSB.m_QueryIdxToTopScoreVec;
		asserta(QueryIdxToTopScoreVec.size() == QSeqCount);
		for (uint qidx = 0; qidx < QSeqCount; ++qidx)
			{
			const string &q = QLabels[qidx];
			const vector<uint16_t> &row = QueryIdxToTopScoreVec[qidx];
			for (uint tidx = 0; tidx < TSeqCount; ++tidx)
				{
				uint16_t score = row[tidx];
				if (score > 0)
					{
					const string &t = TLabels[tidx];
					fprintf(f, "%s\t%s\t%u\n",
						q.c_str(), t.c_str(), score);
					}
				}
			}
		CloseStdioFile(f);
		}
#endif

	if (optset_output3)
		{
		vector<string> QLabels;
		vector<string> TLabels;
		for (uint i = 0; i < QSeqCount; ++i)
			QLabels.push_back(QDB.GetLabel(i));
		for (uint i = 0; i < TSeqCount; ++i)
			TLabels.push_back(TDB.GetLabel(i));
		FILE *fTsv = CreateStdioFile(opt(output3));
		prefilter_kappa::m_RSB.ToLabelsTsv(fTsv, QLabels, TLabels);
		CloseStdioFile(s_fTsv);
		}
	}
