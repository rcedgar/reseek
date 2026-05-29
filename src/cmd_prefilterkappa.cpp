#include "myutils.h"
#include "prefilter_kappa.h"
#include "kappa_mermx.h"
#include "kappa_dex.h"
#include "flat_params.h"
#include "bitdope.h"
#include "lookup.h"

static atomic<uint> s_NextTIdx = 0;
static const kappa_mermx *s_ptrScoreMx;
static const SeqDB *s_ptrQDB = 0;
static const SeqDB *s_ptrTDB = 0;
static const kappa_dex *s_ptrQKmerIndex = 0;
static FILE *s_fTsv = 0;
static time_t s_TimeLastProgress;

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
		uint TSeqIdx = s_NextTIdx++;
		if (TSeqIdx > 0 && TSeqIdx + 1 < TSeqCount)
			{
			time_t now = time(0);
			if (now > s_TimeLastProgress)
				ProgressStep(TSeqIdx, TSeqCount, "Filtering");
			s_TimeLastProgress = now;
			}
		if (TSeqIdx >= TSeqCount)
			return;

		Pref.m_TSeqIdx = TSeqIdx;
		const byte *TSeq = s_ptrTDB->GetByteSeq(TSeqIdx);
		const string &TLabel = s_ptrTDB->GetLabel(TSeqIdx);
		uint TL = s_ptrTDB->GetSeqLength(TSeqIdx);
		Pref.Search(TSeqIdx, TLabel, TSeq, TL);
		}
	}

static void write_tsv(const string &fn)
	{
	if (fn == "") return;

	FILE *fTsv = CreateStdioFile(fn);
	prefilter_kappa::m_RSB.ToTsv(fTsv);
	CloseStdioFile(s_fTsv);
	}

static void write_tsv_with_labels(
	const string &fn,
	const SeqDB &QDB,
	const SeqDB &TDB)
	{
	if (fn == "") return;

	const uint QSeqCount = QDB.GetSeqCount();
	const uint TSeqCount = TDB.GetSeqCount();

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

static void write_tsv_with_scores(
	const string &fn,
	const SeqDB &QDB,
	const SeqDB &TDB)
	{
#if STORE_PAIR_SCORES
	if (fn == "") return;
	const uint QSeqCount = QDB.GetSeqCount();
	const uint TSeqCount = TDB.GetSeqCount();

	vector<string> QLabels;
	vector<string> TLabels;

	for (uint i = 0; i < QSeqCount; ++i)
		QLabels.push_back(QDB.GetLabel(i));
	for (uint i = 0; i < TSeqCount; ++i)
		TLabels.push_back(TDB.GetLabel(i));

	FILE *f = CreateStdioFile(fn);
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
#else
	Die("write_tsv_with_scores() STORE_PAIR_SCORES=0");
#endif
	}

static void bench(
	uint filter_secs,
	const SeqDB &QDB,
	const SeqDB &TDB)
	{
	if (!optset_dope) return;

	const uint QSeqCount = QDB.GetSeqCount();
	const uint TSeqCount = TDB.GetSeqCount();

	bitdope dope;
	lookup look;

	const string &lookupfn =
		(optset_lookup ? opt(lookup) : "../data/scop40x.lookup");
	look.from_tsv(lookupfn);
	dope.m_look = &look;
	dope.from_file(opt(dope));
	dope.set_square();
	ProgressLog("dope %s hits\n", FloatToStr(dope.m_nhit));

	uint npass = 0;
	uint nindope = 0;
	const vector<vector<uint16_t> > &QueryIdxToTopScoreVec =
		prefilter_kappa::m_RSB.m_QueryIdxToTopScoreVec;
	asserta(QueryIdxToTopScoreVec.size() == QSeqCount);
	for (uint qidx = 0; qidx < QSeqCount; ++qidx)
		{
		const string &q = QDB.GetLabel(qidx);
		uint qdomidx = look.get_domidx(q);
		const vector<uint16_t> &row = QueryIdxToTopScoreVec[qidx];
		for (uint tidx = 0; tidx < TSeqCount; ++tidx)
			{
			uint16_t score = row[tidx];
			if (score > 0)
				{
				const string &t = TDB.GetLabel(tidx);
				uint tdomidx = look.get_domidx(t);
				++npass;
				if (dope.in_square_ij(qidx, tidx))
					++nindope;
				}
			}
		}
	double pct = GetPct(nindope, 2*dope.m_nhit);

	Progress("pct=%.1f", pct);
	Progress(" secs=%u", filter_secs);
	Progress(" pattern=%s", flat_params::m_kappa_pattern.c_str());
	Progress(" kmer=%d", flat_params::m_kappa_min_mindiagscore);
	Progress(" diag=%d", flat_params::m_kappa_min_mindiagscore);
	Progress(" npass=%u", npass);
	Progress("\n");

	Log("@FEV@");
	Log("\tpct=%.1f", pct);
	Log("\tsecs=%u", filter_secs);
	Log("\tpattern=%s", flat_params::m_kappa_pattern.c_str());
	Log("\tkmer=%d", flat_params::m_kappa_min_mindiagscore);
	Log("\tdiag=%d", flat_params::m_kappa_min_mindiagscore);
	Log("\tnpass=%u", npass);
	Log("\n");
	}

void cmd_prefilter_kappa()
	{
	const string &QueryKappa_FN = g_Arg1;
	const string &DBFN = opt(db);

	SeqDB QDB;
	SeqDB TDB;

	QDB.FromFasta(QueryKappa_FN);
	TDB.FromFasta(DBFN);

	QDB.ToLetters(g_CharToLetterMu);
	TDB.ToLetters(g_CharToLetterMu);
	const uint QSeqCount = QDB.GetSeqCount();
	const uint TSeqCount = TDB.GetSeqCount();

	void SetQueryNeighborhood(uint QSeqCount);
	SetQueryNeighborhood(QSeqCount);

	prefilter_kappa::init_kappa();
	prefilter_kappa::m_RSB.Init(QSeqCount);

	kappa_dex QKmerIndex;
	QKmerIndex.Init();

	const uint k = flat_params::m_kappa_kmer_nrones;
	const kappa_mermx &GetKappaMerMx(uint k);
	const kappa_mermx &ScoreMx = GetKappaMerMx(k);
	asserta(ScoreMx.m_k == k);

	QKmerIndex.m_KmerSelfScores = ScoreMx.BuildSelfScores_Kmers();
	QKmerIndex.m_MinKmerSelfScore =  flat_params::m_kappa_min_mindiagscore;
	QKmerIndex.FromSeqDB(QDB);
#if DEBUG
	QKmerIndex.Validate();
#endif
	asserta(QKmerIndex.m_k == k);
	asserta(QKmerIndex.m_DictSize == flat_params::m_kappa_dict_size);
	asserta(ScoreMx.m_AS_pow[k] == QKmerIndex.m_DictSize);

	s_ptrQDB = &QDB;
	s_ptrTDB = &TDB;
	s_ptrScoreMx = &ScoreMx;
	s_ptrQKmerIndex = &QKmerIndex;

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

	time_t t_end = time(0);
	uint filter_secs = uint(t_end - t_start);

	uint total = prefilter_kappa::m_RSB.TruncateAllQueryVecs();
	ProgressLog("Prefilter hits  %s\n", FloatToStr(total));

	write_tsv(opt(output));
	write_tsv_with_scores(opt(output2), QDB, TDB);
	write_tsv_with_labels(opt(output3), QDB, TDB);
	bench(filter_secs, QDB, TDB);
	}
