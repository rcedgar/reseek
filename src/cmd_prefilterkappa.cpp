#include "myutils.h"
#include "dssparams.h"
#include "prefilter_kappa.h"
#include "kappa_mermx.h"
#include "kappa_dex.h"
#include "kappa_prefilter_params.h"
#include "bitdope.h"
#include "lookup.h"
#include <chrono>

/////////////////////////////////////////////////////
// Kappa prefilter
// $src/2025-10_reseek_tune [f8b7229]
// reseek v2.9.i86linux64 [89b34d9]
// C:\src\notebooks\2026-04-09_kappa_two_hit_diag_param_manual_explore.txt
// C:\src\2025-10_reseek_tune\bash\test_kappa_prefilter_new_defaults_2026-04-09.bash
// __________________________________________________  Pattern  Kmer   Diag  PFhits  PFTime
// SEPQ0.1=0.296 SEPQ1=0.399 SEPQ10=0.493 Sum3=1.684 |    1111    28    100     9 M   00:11 <<== set these defaults 2026-04-09
// SEPQ0.1=0.290 SEPQ1=0.399 SEPQ10=0.499 Sum3=1.677 | ......Mu bitdope....    12 M   00:33
// 
// Sum3 from hits (v2.7 verysensitive AND kappa filtered)
/////////////////////////////////////////////////////////

//int DSSParams::m_PrefilterMinKappaKmerPairScore = 28;
//int DSSParams::m_PrefilterMinKappaMinDiagScore = 100;
int DSSParams::m_PrefilterMinKappaKmerPairScore = 50;
int DSSParams::m_PrefilterMinKappaMinDiagScore = 150;

uint DSSParams::m_PrefilterKappaKmerNrOnes = 4;
uint DSSParams::m_PrefilterKappaKmerWidth = 4;
uint DSSParams::m_PrefilterKappaDictSize = myipow(32, 4);
string DSSParams::m_PrefilterKappaPattern = "11010001";

//static uint8_t KappaKmerOnesOffsets[] = {0, 1, 3, 7};
//uint8_t *DSSParams::m_PrefilterKappaKmerOnesOffsets =
//	KappaKmerOnesOffsets;
uint8_t *DSSParams::m_PrefilterKappaKmerOnesOffsets = 0;
/////////////////////////////////////////////////////

static uint s_NextTIdx = 0;
static mutex m_NextTIdxLock;
static const kappa_mermx *s_ptrScoreMx;
static const SeqDB *s_ptrQDB = 0;
static const SeqDB *s_ptrTDB = 0;
static const kappa_dex *s_ptrQKmerIndex = 0;
static FILE *s_fTsv = 0;
static time_t s_TimeLastProgress;

static void fill_pattern_offsets(const string &Str, uint8_t *offsets)
	{
	uint n = 0;
	for (uint i = 0; i < SIZE(Str); ++i)
		{
		char c = Str[i];
		asserta(c == '0' || c == '1');
		if (c == '1')
			offsets[n++] = i;
		}
	}

static uint get_nr_pattern_ones(const string &Str)
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
	lookup look;
	bitdope dope;
	if (optset_dope)
		{
		const string &lookupfn =
			(optset_lookup ? opt(lookup) : "../data/scop40x.lookup");
		look.from_tsv(lookupfn);
		dope.m_look = &look;
		dope.from_file(opt(dope));
		dope.set_square();
		ProgressLog("dope %s hits\n", FloatToStr(dope.m_nhit));
		}

	asserta(optset_logodds);
	const double scalef = (optset_scalef ? opt(scalef) : 10);
	void load_kappa_integer_logodds(const string &fn, double scalef);
	load_kappa_integer_logodds(opt(logodds), scalef);

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

	if (optset_kappa_pattern)
		DSSParams::m_PrefilterKappaPattern = opt(kappa_pattern);
	uint k = get_nr_pattern_ones(DSSParams::m_PrefilterKappaPattern);
	uint K = uint(DSSParams::m_PrefilterKappaPattern.size());
	DSSParams::m_PrefilterKappaKmerOnesOffsets = myalloc(uint8_t, k);
	fill_pattern_offsets(DSSParams::m_PrefilterKappaPattern,
		DSSParams::m_PrefilterKappaKmerOnesOffsets);

	DSSParams::m_PrefilterKappaKmerNrOnes = k; 
	DSSParams::m_PrefilterKappaKmerWidth = K;
	DSSParams::m_PrefilterKappaDictSize = myipow(32, k);

	if (optset_kappa_minkmerscore)
		DSSParams::m_PrefilterMinKappaKmerPairScore = opt(kappa_minkmerscore);
	if (optset_kappa_mindiagscore)
		DSSParams::m_PrefilterMinKappaMinDiagScore = opt(kappa_mindiagscore);

	kappa_dex QKmerIndex;
	QKmerIndex.Init();

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
	uint total = prefilter_kappa::m_RSB.TruncateAllQueryVecs();
	ProgressLog("Prefilter hits  %s\n", FloatToStr(total));

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

	if (optset_dope)
		{
		uint nhit = 0;
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
					++nhit;
					if (dope.in_square_ij(qidx, tidx))
						++nindope;
					}
				}
			}
		double pct = GetPct(nindope, 2*dope.m_nhit);

		//ProgressLog("%u / %u filter hits also in dope\n",
		//	nindope, nhit);
		//ProgressLog("%u / %u dope passed filter (%.2f%%)\n",
		//	nindope, 2*dope.m_nhit, pct);

		Progress("pct=%.1f", pct);
		Progress(" secs=%u", filter_secs);
		Progress(" pattern=%s", DSSParams::m_PrefilterKappaPattern.c_str());
		Progress(" kmer=%d", DSSParams::m_PrefilterMinKappaKmerPairScore);
		Progress(" diag=%d", DSSParams::m_PrefilterMinKappaMinDiagScore);
		Progress("\n");

		Log("@FEV@");
		Log("\tpct=%.1f", pct);
		Log("\tsecs=%u", filter_secs);
		Log("\tpattern=%s", DSSParams::m_PrefilterKappaPattern.c_str());
		Log("\tkmer=%d", DSSParams::m_PrefilterMinKappaKmerPairScore);
		Log("\tdiag=%d", DSSParams::m_PrefilterMinKappaMinDiagScore);
		Log("\n");
		}
	}
