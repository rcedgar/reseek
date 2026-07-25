#include "myutils.h"
#include "kappa_filter.h"
#include "kappa_hsp.h"
#include "flat_params.h"
#include "seqinfo.h"
#include <mutex>

RankedScoresBag kappa_filter::m_RSB;
uint8_t **kappa_filter::m_query_kappa_codeseq_vec = 0;
const uint *kappa_filter::m_query_lengths = 0;
const vector<string> *kappa_filter::m_query_labels = 0;
kappa_seqsource *kappa_filter::m_db_seqsource = 0;
uint kappa_filter::m_QSeqCount = 0;
atomic<time_t> kappa_filter::m_time_last_progress;

const kappa_mermx *kappa_filter::m_ptrScoreMx;
const kappa_dex *kappa_filter::m_ptrQKmerIndex;
FILE *kappa_filter::m_f_prehsp_dump = 0;
mutex kappa_filter::m_prehsp_dump_mutex;
atomic<uint64> kappa_filter::m_n_targets{0};
atomic<uint64> kappa_filter::m_n_targets_skipped{0};
atomic<uint64> kappa_filter::m_n_index_postings{0};
atomic<uint64> kappa_filter::m_n_prehsp{0};
atomic<uint64> kappa_filter::m_n_hsp{0};
atomic<uint64> kappa_filter::m_n_hsp_prune{0};
atomic<uint64> kappa_filter::m_n_hsp_cells{0};
atomic<uint64> kappa_filter::m_n_rsb_qt_pairs{0};
vector<KappaThreadDiag> kappa_filter::m_thread_diags;
bool g_QueryNeighborhood = true;

void fill_pattern_offsets(const string &Str, uint8_t *offsets)
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

uint get_nr_pattern_ones(const string &Str)
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

bool kappa_filter::m_init_kappa_done = false;
void kappa_filter::init_kappa()
	{
	asserta(!m_init_kappa_done);

	flat_params::init_kappa();
	kappa_filter::m_RSB.m_B = flat_params::m_rsb_size;

	const char *diag_mode = "unique_fine";
	if (flat_params::m_kappa_onehitdiag)
		diag_mode = "onehit_insert";
	else if (flat_params::m_kappa_twohitdiag)
		diag_mode = "twohit_dupes";
	ProgressLog("Kappa diag mode %s\n", diag_mode);

	m_init_kappa_done = true;
	}

//////////////////////////////////////////////
// 	FindHSP searches for the highest-scoring
// 	ungapped alignment on a given diagonal.
//////////////////////////////////////////////
int kappa_filter::FindHSP(const byte *QSeq, uint QL, int Diag) const
	{
	return kappa_find_hsp(QSeq, m_TSeq, int(QL), int(m_TL), Diag);
	}

int kappa_filter::FindHSP(uint QSeqIdx, int Diag) const
	{
	//const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	//const uint QL = m_QDB->GetSeqLength(QSeqIdx);
	const byte *QSeq = m_query_kappa_codeseq_vec[QSeqIdx];
	const uint QL = m_query_lengths[QSeqIdx];
#if TRACE
	if (DoTrace(QSeqIdx))
		Log("FindHSP(QL=%u, Diag=%d)\n", QL, Diag);
#endif
	return FindHSP(QSeq, QL, Diag);
	}

//////////////////////////////////////////////
// 	FindHSP plus "traceback", i.e. returns
// 	start position and length of HSP.
//////////////////////////////////////////////
int kappa_filter::FindHSP2(const byte *QSeq, uint QL, int Diag,
							  int &Lo, int &Len) const
	{
	asserta(Diag >= 0);
	const int LQ = int(QL);
	const int LT = int(m_TL);
	const int d = Diag;
	int mini = LQ - d - 1;
	if (mini < 0)
		mini = 0;
	int minj = d + 1 - LQ;
	if (minj < 0)
		minj = 0;
	int maxi = LQ + LT - d - 2;
	if (maxi >= LQ)
		maxi = LQ - 1;
	const int n = maxi - mini + 1;
	asserta(n > 0);

	const byte *q = QSeq + mini;
	const byte *t = m_TSeq + minj;
	int B = 0;
	int F = 0;
	int CurrLen = 0;
	Lo = 0;
	Len = 0;
	int SuffixLo = 0;
	for (int k = 0; k < n; ++k)
		{
		const byte bq = *q++;
		const byte bt = *t++;
#if !defined(NDEBUG)
		assert(bq < KAPPA_AS);
		assert(bt < KAPPA_AS);
#endif
		const int Score = int(kappa32_flat_logodds[(unsigned) bq*32u + (unsigned) bt]);
		F += Score;
		if (F > B)
			{
			B = F;
			Lo = SuffixLo;
			Len = ++CurrLen;
			}
		else if (F > 0)
			++CurrLen;
		else
			{
			F = 0;
			SuffixLo = k+1;
			CurrLen = 0;
			}
		}
	return B;
	}

int kappa_filter::FindHSP2(uint QSeqIdx, int Diag, int &Lo, int &Len) const
	{
	//const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	//const uint QL = m_QDB->GetSeqLength(QSeqIdx);
	const byte *QSeq = m_query_kappa_codeseq_vec[QSeqIdx];
	const uint QL = m_query_lengths[QSeqIdx];
	return FindHSP2(QSeq, QL, Diag, Lo, Len);
	}

kappa_filter::~kappa_filter()
	{
	if (!m_RSBPending.empty())
		m_RSB.AddScoresBatch(m_RSBPending);
	flush_filter_stats();
	}

void kappa_filter::flush_filter_stats()
	{
	if (m_local_n_targets != 0)
		m_n_targets.fetch_add(m_local_n_targets, memory_order_relaxed);
	if (m_local_n_targets_skipped != 0)
		m_n_targets_skipped.fetch_add(m_local_n_targets_skipped, memory_order_relaxed);
	if (m_local_n_index_postings != 0)
		m_n_index_postings.fetch_add(m_local_n_index_postings, memory_order_relaxed);
	if (m_local_n_prehsp != 0)
		m_n_prehsp.fetch_add(m_local_n_prehsp, memory_order_relaxed);
	if (m_local_n_hsp != 0)
		m_n_hsp.fetch_add(m_local_n_hsp, memory_order_relaxed);
	if (m_local_n_hsp_prune != 0)
		m_n_hsp_prune.fetch_add(m_local_n_hsp_prune, memory_order_relaxed);
	if (m_local_n_hsp_cells != 0)
		m_n_hsp_cells.fetch_add(m_local_n_hsp_cells, memory_order_relaxed);
	if (m_local_n_rsb_qt_pairs != 0)
		m_n_rsb_qt_pairs.fetch_add(m_local_n_rsb_qt_pairs, memory_order_relaxed);
	m_local_n_targets = 0;
	m_local_n_targets_skipped = 0;
	m_local_n_index_postings = 0;
	m_local_n_prehsp = 0;
	m_local_n_hsp = 0;
	m_local_n_hsp_prune = 0;
	m_local_n_hsp_cells = 0;
	m_local_n_rsb_qt_pairs = 0;
	}

void kappa_filter::record_thread_diag(uint threadidx, uint batches, uint secs)
	{
	asserta(threadidx < SIZE(m_thread_diags));
	KappaThreadDiag &d = m_thread_diags[threadidx];
	d.batches = batches;
	d.targets = m_local_n_targets;
	d.hsp = m_local_n_hsp;
	d.hsp_cells = m_local_n_hsp_cells;
	d.ticks_seed = m_local_ticks_seed;
	d.ticks_diag = m_local_ticks_diag;
	d.ticks_hsp = m_local_ticks_hsp;
	d.ticks_rsb = m_local_ticks_rsb;
	d.secs = secs;
	}

void kappa_filter::reset_filter_stats()
	{
	m_n_targets.store(0, memory_order_relaxed);
	m_n_targets_skipped.store(0, memory_order_relaxed);
	m_n_index_postings.store(0, memory_order_relaxed);
	m_n_prehsp.store(0, memory_order_relaxed);
	m_n_hsp.store(0, memory_order_relaxed);
	m_n_hsp_prune.store(0, memory_order_relaxed);
	m_n_hsp_cells.store(0, memory_order_relaxed);
	m_n_rsb_qt_pairs.store(0, memory_order_relaxed);
	}

void kappa_filter::alloc()
	{
	asserta(m_QSeqCount > 0);
	m_RSBPending.clear();
	m_RSBPending.reserve(RSB_BATCH);

	m_QSeqIdxToBestDiagScore = myalloc(uint16_t, m_QSeqCount);
	m_QSeqIdxsWithTwoHitDiag = myalloc(uint16_t, m_QSeqCount);

	for (uint i = 0; i < m_QSeqCount; ++i)
		{
		m_QSeqIdxToBestDiagScore[i] = 0;
		m_QSeqIdxsWithTwoHitDiag[i] = UINT16_MAX;
		}

	bool TargetNeighborhood = !g_QueryNeighborhood;
	if (TargetNeighborhood)
		m_NeighborKmers = myalloc(uint, flat_params::m_kappa_dict_size);
	else
		m_NeighborKmers = 0;
	m_NrQueriesWithTwoHitDiag = 0;
	}

void kappa_filter::SetQDB(const SeqDB &QDB)
	{
	Die("kappa_filter::SetQDB()");
	//asserta(m_init_kappa_done);

	//m_QDB = &QDB;
	//m_QSeqCount = QDB.GetSeqCount();

	//m_RSBPending.clear();
	//m_RSBPending.reserve(RSB_BATCH);

	//m_QSeqIdxToBestDiagScore = myalloc(uint16_t, m_QSeqCount);
	//m_QSeqIdxsWithTwoHitDiag = myalloc(uint16_t, m_QSeqCount);

	//for (uint i = 0; i < m_QSeqCount; ++i)
	//	{
	//	m_QSeqIdxToBestDiagScore[i] = 0;
	//	m_QSeqIdxsWithTwoHitDiag[i] = UINT16_MAX;
	//	}

	//bool TargetNeighborhood = !g_QueryNeighborhood;
	//if (TargetNeighborhood)
	//	m_NeighborKmers = myalloc(uint, flat_params::m_kappa_dict_size);
	//else
	//	m_NeighborKmers = 0;
	//m_NrQueriesWithTwoHitDiag = 0;
	}

void kappa_filter::Search_TargetKmers()
	{
	m_QKmerIndex->GetKmers(m_TSeq, m_TL, m_TKmers);
	const uint NK = SIZE(m_TKmers);
	if (g_QueryNeighborhood)
		{
		for (uint TPos = 0; TPos < NK; ++TPos)
			{
			uint TKmer = m_TKmers[TPos];
			if (TKmer != UINT_MAX)
				{
#if TRACE
				m_TBaseKmer = TKmer;
#endif
				Search_TargetKmer(TKmer, TPos);
				}
			}
		}
	else
		{
		for (uint TPos = 0; TPos < NK; ++TPos)
			{
			uint Kmer = m_TKmers[TPos];
			if (Kmer != UINT_MAX)
				Search_TargetKmerNeighborhood(Kmer, TPos);
			}
		}
	}

void kappa_filter::Search_TargetSeq(uint TSeqIdx, const string &TLabel,
				   const byte *TSeq, uint TL)
	{
	m_TSeqIdx = TSeqIdx;
	m_TLabel = TLabel;
	m_TSeq = TSeq;
	m_TL = TL;

	Reset();

	TICKS t0 = GetClockTicks();
	Search_TargetKmers();
	TICKS t1 = GetClockTicks();
	m_local_ticks_seed += t1 - t0;

	if (flat_params::m_kappa_onehitdiag)
		{
		DumpPreHSPHits();
		t0 = GetClockTicks();
		ExtendOneHitDiagsToHSPs();
		t1 = GetClockTicks();
		m_local_ticks_hsp += t1 - t0;
		}
	else
		{
		t0 = GetClockTicks();
		FindTwoHitDiags();
		t1 = GetClockTicks();
		m_local_ticks_diag += t1 - t0;

		DumpPreHSPHits();

		t0 = GetClockTicks();
		ExtendTwoHitDiagsToHSPs();
		t1 = GetClockTicks();
		m_local_ticks_hsp += t1 - t0;
		}
	}

void kappa_filter::Search_TargetKmerNeighborhood(uint Kmer, uint TPos)
	{
	if (Kmer == UINT_MAX)
		return;
#if TRACE
	m_TBaseKmer = Kmer;
#endif
	assert(Kmer < flat_params::m_kappa_dict_size);
	if (m_KmerSelfScores[Kmer] < flat_params::m_kappa_min_kmerpairscore)
		return;
	short MinKmerScore =  flat_params::m_kappa_min_kmerpairscore;

// Construct high-scoring neighborhood
	const uint HSKmerCount =
		m_ScoreMx->GetHighScoringKmers(Kmer, MinKmerScore, m_NeighborKmers);

#if TRACE
	string Tmp;
	Log("Search_TargetKmerNeighborhood TPos=%u Kmer=%s minscore=%d |nbrs|=%d\n",
		TPos, KmerToStr(Kmer, Tmp), MinKmerScore, HSKmerCount);
#endif
	for (uint HSKmerIdx = 0; HSKmerIdx < HSKmerCount; ++HSKmerIdx)
		{
		uint HSKmer = m_NeighborKmers[HSKmerIdx];
		Search_TargetKmer(HSKmer, TPos);
		}
	}

void kappa_filter::Search_TargetKmer(uint TKmer, uint TPos)
	{
	uint64_t RowSize = m_QKmerIndex->GetRowSize(TKmer);
#if TRACE
	{
	string KmerStr;
	m_QKmerIndex->KmerToStr(m_TBaseKmer, KmerStr);
	Log("Search_TargetKmer(TPos=%u, TKmer=%s) RowSize=%s\n",
		TPos, KmerStr.c_str(), Int64ToStr(RowSize));
	}
#endif
	if (RowSize == 0)
		return;
	m_local_n_index_postings += RowSize;
	uint64_t DataOffset = m_QKmerIndex->GetRowStart(TKmer);
	for (uint64_t ColIdx = 0; ColIdx < RowSize; ++ColIdx)
		{
		uint32_t QSeqIdx;
		uint16_t QSeqPos;
		m_QKmerIndex->Get(DataOffset++, QSeqIdx, QSeqPos);
		asserta(QSeqIdx < m_QSeqCount);
		//uint QL32 = m_QDB->GetSeqLength(QSeqIdx);
		uint QL32 = m_query_lengths[QSeqIdx];
		asserta(QL32 < UINT16_MAX);
		uint16_t Diag = uint16_t(QL32 + TPos - QSeqPos - 1);
#if TRACE
		{
		string TKmerStr;
		string QKmerStr;
		uint QKmer = GetQKmer(QSeqIdx, QSeqPos);
		m_QKmerIndex->KmerToStr(m_TBaseKmer, TKmerStr);
		m_QKmerIndex->KmerToStr(QKmer, QKmerStr);
		const MerMx &MM = GetKappaMerMx(m_QKmerIndex->m_k);
		int KmerPairScore = MM.GetScoreKmerPair(TKmer, QKmer);

		Log("@K@  [%4u] %5s  [%4u] %5s  /%5u/  %+3d\n",
			QSeqPos, QKmerStr.c_str(),
			TPos, TKmerStr.c_str(),
			Diag, KmerPairScore);
		}
#endif
		if (Diag > m_Mask14)
			continue;
		if (flat_params::m_kappa_onehitdiag)
			OneHitDiagAdd(QSeqIdx, Diag);
		else
			m_DiagBag.Add(QSeqIdx, Diag);
		}
	}
	
void kappa_filter::FindTwoHitDiags()
	{
	if (flat_params::m_kappa_twohitdiag)
		m_DiagBag.SetDupes();
	else
		m_DiagBag.SetUniqueFine();
#if DEBUG
	//m_DiagBag.Validate(m_QSeqCount, INT16_MAX);
#endif
	}

void kappa_filter::OneHitDiagAdd(uint SeqIdx, uint16_t Diag)
	{
	asserta(SeqIdx < UINT16_MAX);
	const uint32_t pair = (uint32_t(SeqIdx) << 16) | uint32_t(Diag);
	m_OneHitDiags.insert(pair);
	}

void kappa_filter::ExtendOneHitDiagsToHSPs()
	{
	m_NrQueriesWithTwoHitDiag = 0;
	m_local_n_prehsp += m_OneHitDiags.size();
	for (set<uint32_t>::const_iterator iter = m_OneHitDiags.begin();
		 iter != m_OneHitDiags.end(); ++iter)
		{
		const uint32_t pair = *iter;
		const uint32_t QSeqIdx = (pair >> 16);
		const uint16_t Diag = uint16_t(pair & 0xffff);
		const int DiagScore = ExtendDiagToHSP(QSeqIdx, Diag);
		AddTwoHitDiag(QSeqIdx, Diag, DiagScore);
		}
	}

void kappa_filter::set_prehsp_dump(FILE *f, const vector<string> *query_labels)
	{
	m_f_prehsp_dump = f;
	m_query_labels = query_labels;
	}

void kappa_filter::write_prehsp_tsv_header(FILE *f, const char *path_tag)
	{
	asserta(f != 0);
	fprintf(f, "# idx_prefilter_prehsp\n");
	fprintf(f, "# path\t%s\n", path_tag);
	const char *diag_mode = "unique_fine";
	if (flat_params::m_kappa_onehitdiag)
		diag_mode = "onehit_insert";
	else if (flat_params::m_kappa_twohitdiag)
		diag_mode = "twohit_dupes";
	fprintf(f, "# diag_mode\t%s\n", diag_mode);
	fprintf(f, "# minkmerscore\t%d\n", flat_params::m_kappa_min_kmerpairscore);
	fprintf(f, "# fields\tq_label\tt_label\tqidx\ttidx\tdiag\n");
	}

void kappa_filter::write_prehsp_hit(FILE *f,
	const char *qlabel, const char *tlabel,
	uint qidx, uint tidx, uint16_t diag)
	{
	asserta(f != 0);
	if (qlabel == 0) qlabel = "-";
	if (tlabel == 0) tlabel = "-";
	fprintf(f, "%s\t%s\t%u\t%u\t%u\n", qlabel, tlabel, qidx, tidx, diag);
	}

void kappa_filter::DumpPreHSPHits() const
	{
	if (m_f_prehsp_dump == 0)
		return;

	const char *tlabel = m_TLabel.c_str();
	lock_guard<mutex> lock(m_prehsp_dump_mutex);

	if (flat_params::m_kappa_onehitdiag)
		{
		for (set<uint32_t>::const_iterator iter = m_OneHitDiags.begin();
			 iter != m_OneHitDiags.end(); ++iter)
			{
			const uint32_t pair = *iter;
			const uint qidx = (pair >> 16);
			const uint16_t diag = uint16_t(pair & 0xffff);
			const char *qlabel = "-";
			if (m_query_labels != 0 && qidx < SIZE(*m_query_labels))
				qlabel = (*m_query_labels)[qidx].c_str();
			write_prehsp_hit(m_f_prehsp_dump, qlabel, tlabel,
				qidx, m_TSeqIdx, diag);
			}
		return;
		}

	const uint n = m_DiagBag.m_DupeCount;
	for (uint i = 0; i < n; ++i)
		{
		const uint qidx = m_DiagBag.m_DupeSeqIdxs[i];
		const uint16_t diag = m_DiagBag.m_DupeDiags[i];
		const char *qlabel = "-";
		if (m_query_labels != 0 && qidx < SIZE(*m_query_labels))
			qlabel = (*m_query_labels)[qidx].c_str();
		write_prehsp_hit(m_f_prehsp_dump, qlabel, tlabel,
			qidx, m_TSeqIdx, diag);
		}
	}

void kappa_filter::GetResults(vector<uint> &QSeqIdxs,
						   vector<uint16_t> &DiagScores) const
	{
	QSeqIdxs.clear();
	DiagScores.clear();
	DiagScores.reserve(m_NrQueriesWithTwoHitDiag);
 	for (uint i = 0; i < m_NrQueriesWithTwoHitDiag; ++i)
		{
		uint QSeqIdx = m_QSeqIdxsWithTwoHitDiag[i];
		uint16_t DiagScore = m_QSeqIdxToBestDiagScore[QSeqIdx];

		QSeqIdxs.push_back(QSeqIdx);
		DiagScores.push_back(DiagScore);
		}
	}

void kappa_filter::AddTwoHitDiag(uint QSeqIdx, uint16_t Diag, int DiagScore)
	{
	if (DiagScore <= 0)
		return;
	if (DiagScore < flat_params::m_kappa_min_diagscore)
		return;
	asserta(QSeqIdx < UINT16_MAX);
	if (DiagScore >= UINT16_MAX)
		DiagScore = UINT16_MAX-1;
	uint16_t BestDiagScoreT = m_QSeqIdxToBestDiagScore[QSeqIdx];
	if (BestDiagScoreT == 0)
		{
		m_QSeqIdxsWithTwoHitDiag[m_NrQueriesWithTwoHitDiag++] = QSeqIdx;
		m_QSeqIdxToBestDiagScore[QSeqIdx] = DiagScore;
#if TRACE
		if (DoTrace(QSeqIdx)) Log("AddTwoHitDiag(TSeqIdx=%u, QSeqIdx=%u, Diag=%u, DiagScore=%d) (first)\n",
								  m_TSeqIdx, QSeqIdx, Diag, DiagScore);
#endif
		}
	else if (DiagScore > m_QSeqIdxToBestDiagScore[QSeqIdx])
		{
		m_QSeqIdxToBestDiagScore[QSeqIdx] = DiagScore;
#if TRACE
		if (DoTrace(QSeqIdx)) Log("AddTwoHitDiag(TSeqIdx=%u, QSeqIdx=%u, Diag=%u, DiagScore=%d) (better)\n",
								  m_TSeqIdx, QSeqIdx, Diag, DiagScore);
#endif
		}
	}

void kappa_filter::ExtendTwoHitDiagsToHSPs()
	{
	const uint DupeCount = m_DiagBag.m_DupeCount;
	m_local_n_prehsp += DupeCount;
	m_NrQueriesWithTwoHitDiag = 0;
	for (uint i = 0; i < DupeCount; ++i)
		{
		const uint32_t QSeqIdx = m_DiagBag.m_DupeSeqIdxs[i];
		const uint16_t Diag = m_DiagBag.m_DupeDiags[i];
		const int DiagScore = ExtendDiagToHSP(QSeqIdx, Diag);
		AddTwoHitDiag(QSeqIdx, Diag, DiagScore);
		}
	}

int kappa_filter::ExtendDiagToHSP(uint32_t QSeqIdx, uint16_t Diag)
	{
	const byte *QSeq = m_query_kappa_codeseq_vec[QSeqIdx];
	const uint QL = m_query_lengths[QSeqIdx];

	int mini, minj, n;
	kappa_get_hsp_limits(int(QL), int(m_TL), int(Diag),
		mini, minj, n);

	if (flat_params::m_kappa_hsp_rsb_prune &&
		m_RSB.m_AnyLoScoreActive.load(std::memory_order_relaxed))
		{
		const uint16_t LoScore = m_RSB.GetLoScore(QSeqIdx);
		if (LoScore > 0)
			{
			const int bound = n * flat_params::m_kappa_max_pos_logodds;
			if (bound < LoScore)
				{
				++m_local_n_hsp_prune;
				return 0;
				}
			}
		}

	++m_local_n_hsp;
	m_local_n_hsp_cells += uint64(n);
	int DiagScore = FindHSP(QSeq, QL, Diag);
#if TRACE
	LogDiag(QSeqIdx, Diag);
#endif
	return DiagScore;
	}

void kappa_filter::Reset()
	{
	for (uint HitIdx = 0; HitIdx < m_NrQueriesWithTwoHitDiag; ++HitIdx)
		{
		uint QSeqIdx = m_QSeqIdxsWithTwoHitDiag[HitIdx];
		m_QSeqIdxToBestDiagScore[QSeqIdx] = 0;
		}
#if DEBUG
	{
	for (uint SeqIdx = 0; SeqIdx < m_QSeqCount; ++SeqIdx)
		{
		assert(m_QSeqIdxToBestDiagScore[SeqIdx] == 0);
		}
	}
#endif
	m_NrQueriesWithTwoHitDiag = 0;
	m_DiagBag.Reset();
	m_OneHitDiags.clear();
	}

void kappa_filter::LogDiag(uint QSeqIdx, uint16_t Diag) const
	{
	//const byte *QSeq = m_QDB->GetByteSeq(QSeqIdx);
	//uint QL = m_QDB->GetSeqLength(QSeqIdx);
	const byte *QSeq = m_query_kappa_codeseq_vec[QSeqIdx];
	const uint QL = m_query_lengths[QSeqIdx];
	string QSeq_ascii;
	for (uint i = 0; i < QL; ++i)
		QSeq_ascii += g_LetterToCharMu[QSeq[i]];
	int Score = FindHSP(QSeq, QL, Diag);
	int Lo, Len;
	int Score2 = FindHSP2(QSeq, QL, Diag, Lo, Len);
	//const string &QLabel = m_QDB->GetLabel(QSeqIdx);
	Log("LogDiag(%u) lo %d, len %d, score %d\n",
		Diag, Lo, Len, Score);
	diag dg(QL, m_TL);
	int ilo = dg.getmini(Diag) + Lo;
	int jlo = dg.getminj(Diag) + Lo;
	string TSeq_ascii;
	for (uint i = 0; i < m_TL; ++i)
		TSeq_ascii += g_LetterToCharMu[m_TSeq[i]];
	Log(" Q %*.*s\n", Len, Len, QSeq_ascii.c_str() + ilo);
	Log(" T %*.*s\n", Len, Len, TSeq_ascii.c_str() + jlo);
	asserta(Score2 == Score);
	}

void kappa_filter::Search(uint TSeqIdx, const string &TLabel,
				const byte *TSeq, uint TL)
	{
	++m_local_n_targets;
	Search_TargetSeq(TSeqIdx, TLabel, TSeq, TL);

	TICKS t0 = GetClockTicks();
	m_local_n_rsb_qt_pairs += m_NrQueriesWithTwoHitDiag;
	for (uint i = 0; i < m_NrQueriesWithTwoHitDiag; ++i)
		{
		uint QSeqIdx = m_QSeqIdxsWithTwoHitDiag[i];
		uint16_t DiagScore = m_QSeqIdxToBestDiagScore[QSeqIdx];
		RankedScoreBatchEntry e;
		e.QueryIdx = QSeqIdx;
		e.TargetIdx = m_TSeqIdx;
		e.Score = DiagScore;
		m_RSBPending.push_back(e);
		if (m_RSBPending.size() >= RSB_BATCH)
			m_RSB.AddScoresBatch(m_RSBPending);
		}
	TICKS t1 = GetClockTicks();
	m_local_ticks_rsb += t1 - t0;
	}

uint kappa_filter::GetQKmer(uint QSeqIdx, uint QPos) const
	{
	//const byte *Q = m_QDB->GetByteSeq(QSeqIdx);
	const byte *Q = m_query_kappa_codeseq_vec[QSeqIdx];
	uint Kmer = m_QKmerIndex->BytesToKmer(Q + QPos);
	return Kmer;
	}

void kappa_filter::LogQueryKmers(uint QSeqIdx) const
	{
	//const byte *Q = m_QDB->GetByteSeq(QSeqIdx);
	//const uint QL = m_QDB->GetSeqLength(QSeqIdx);
	const byte *QSeq = m_query_kappa_codeseq_vec[QSeqIdx];
	const uint QL = m_query_lengths[QSeqIdx];
	Log("\n");
	Log("kappa_filter::LogQueryKmers() QL=%u\n", QL);
	for (uint PosQ = 0; PosQ + kappa_dex::m_K <= QL; ++PosQ)
		{
		uint Kmer = m_QKmerIndex->BytesToKmer(QSeq + PosQ);
		string tmp;
		const char *KmerStr = m_QKmerIndex->KmerToStr(Kmer, tmp);
		Log("[%4u]  %08x  %s\n", PosQ, Kmer, KmerStr);
		}
	}

void kappa_filter::LogTargetKmers() const
	{
	Log("\n");
	Log("kappa_filter::LogTargetKmers() TL=%u >%s\n", 
		m_TL, m_TLabel);
	for (uint PosT = 0; PosT + kappa_dex::m_K <= m_TL; ++PosT)
		{
		uint Kmer = m_QKmerIndex->BytesToKmer(m_TSeq + PosT);
		string tmp;
		const char *KmerStr = m_QKmerIndex->KmerToStr(Kmer, tmp);
		Log("[%4u]  %08x  %s\n", PosT, Kmer, KmerStr);
		}
	}

void kappa_filter::static_thread_body(uint threadidx)
	{
	asserta(kappa_filter::m_QSeqCount > 0);

	const time_t t_thread_start = time(0);

	kappa_filter Pref;
	Pref.m_ScoreMx = m_ptrScoreMx;
	Pref.m_QKmerIndex = m_ptrQKmerIndex;
	Pref.m_KmerSelfScores = m_ptrQKmerIndex->m_KmerSelfScores;
	Pref.alloc();

	ObjMgr OM;

	uint counter = 0;
	for (;;)
		{
		SeqInfo *TargetSI = OM.GetSeqInfo();
		bool ok = m_db_seqsource->GetNext(TargetSI);
		if (!ok)
			break;
		if ((counter++)%10 == 0)
			{
			time_t now = time(0);
			if (now > m_time_last_progress)
				{
				static mutex s_progress_lock;
				s_progress_lock.lock();
				uint pctx10 = m_db_seqsource->GetPctDoneX10();
				if (pctx10 >= 999) pctx10 = 998;
				ProgressStep(pctx10, 1000, "Kappa filter");
				s_progress_lock.unlock();
				m_time_last_progress = now;
				}
			}

		uint TL = TargetSI->m_L;
		if (TL < flat_params::m_kappa_min_chainlength)
			{
			++Pref.m_local_n_targets_skipped;
			OM.Down(TargetSI);
			continue;
			}

		const uint TSeqIdx = TargetSI->m_Index;
		const uint8_t * TSeq = TargetSI->m_Seq;
		const string &TLabel = TargetSI->m_Label;
		Pref.Search(TSeqIdx, TLabel, TSeq, TL);

		OM.Down(TargetSI);
		}

	const uint secs = uint(time(0) - t_thread_start);
	Pref.record_thread_diag(threadidx, 0, secs);
	}

void kappa_filter::static_bcb_thread_body(uint threadidx)
	{
	asserta(kappa_filter::m_QSeqCount > 0);
	asserta(m_db_seqsource != 0);
	asserta(m_db_seqsource->m_KSSS == KSSS_bcb);

	const time_t t_thread_start = time(0);
	kappa_seqsource &db = *m_db_seqsource;

	kappa_filter Pref;
	Pref.m_ScoreMx = m_ptrScoreMx;
	Pref.m_QKmerIndex = m_ptrQKmerIndex;
	Pref.m_KmerSelfScores = m_ptrQKmerIndex->m_KmerSelfScores;
	Pref.alloc();

	uint counter = 0;
	uint batches = 0;
	for (;;)
		{
		KssBcbBatch *batch = db.claim_bcb_batch();
		if (batch == 0)
			break;
		++batches;

		for (uint i = 0; i < batch->count; ++i)
			{
			const KssBcbSlot &slot = batch->slots[i];
			asserta(slot.label != 0);

			++db.m_bcb_done_count;

			if ((counter++)%10 == 0)
				{
				time_t now = time(0);
				if (now > m_time_last_progress)
					{
					static mutex s_progress_lock;
					s_progress_lock.lock();
					uint pctx10 = db.GetPctDoneX10();
					if (pctx10 >= 999) pctx10 = 998;
					ProgressStep(pctx10, 1000, "Kappa filter");
					s_progress_lock.unlock();
					m_time_last_progress = now;
					}
				}

			if (slot.L < flat_params::m_kappa_min_chainlength)
				{
				++Pref.m_local_n_targets_skipped;
				continue;
				}

			Pref.Search(slot.idx, *slot.label,
				slot.kappa.data(), slot.L);
			}

		db.release_bcb_batch(batch);
		}

	const uint secs = uint(time(0) - t_thread_start);
	Pref.record_thread_diag(threadidx, batches, secs);
	}

void kappa_filter::log_filter_diagnostics(uint NQ, uint NDB, uint total, uint rsb_sat)
	{
	const uint B = kappa_filter::m_RSB.m_B;
	uint64 idx_postings = 0;
	if (m_ptrQKmerIndex != 0)
		idx_postings = m_ptrQKmerIndex->m_Size;

	ProgressLog("Kappa prefilter hits  %s\n", FloatToStr(total));
	ProgressLog("Kappa filter stats  NQ=%u NDB=%u B=%u idx_postings=%llu\n",
		NQ, NDB, B, (unsigned long long) idx_postings);
	ProgressLog("  targets=%llu skipped=%llu postings_visited=%llu prehsp=%llu\n",
		(unsigned long long) m_n_targets.load(memory_order_relaxed),
		(unsigned long long) m_n_targets_skipped.load(memory_order_relaxed),
		(unsigned long long) m_n_index_postings.load(memory_order_relaxed),
		(unsigned long long) m_n_prehsp.load(memory_order_relaxed));
	ProgressLog("  hsp=%llu hsp_prune=%llu hsp_cells=%llu rsb_qt=%llu rsb_sat=%u\n",
		(unsigned long long) m_n_hsp.load(memory_order_relaxed),
		(unsigned long long) m_n_hsp_prune.load(memory_order_relaxed),
		(unsigned long long) m_n_hsp_cells.load(memory_order_relaxed),
		(unsigned long long) m_n_rsb_qt_pairs.load(memory_order_relaxed),
		rsb_sat);

	const uint nthreads = SIZE(m_thread_diags);
	uint busy = 0;
	uint batches_total = 0;
	uint max_secs = 0;
	uint min_busy_secs = UINT_MAX;
	uint64 sum_secs = 0;
	uint64 ticks_seed = 0;
	uint64 ticks_diag = 0;
	uint64 ticks_hsp = 0;
	uint64 ticks_rsb = 0;
	for (uint ti = 0; ti < nthreads; ++ti)
		{
		const KappaThreadDiag &d = m_thread_diags[ti];
		batches_total += d.batches;
		ticks_seed += d.ticks_seed;
		ticks_diag += d.ticks_diag;
		ticks_hsp += d.ticks_hsp;
		ticks_rsb += d.ticks_rsb;
		if (d.targets == 0 && d.batches == 0)
			continue;
		++busy;
		sum_secs += d.secs;
		if (d.secs > max_secs)
			max_secs = d.secs;
		if (d.secs < min_busy_secs)
			min_busy_secs = d.secs;
		}
	if (busy == 0)
		min_busy_secs = 0;

	ProgressLog("Kappa thread balance  nthreads=%u busy=%u batches=%u max_secs=%u min_busy_secs=%u sum_secs=%llu\n",
		nthreads, busy, batches_total, max_secs, min_busy_secs,
		(unsigned long long) sum_secs);

	const uint64 ticks_all = ticks_seed + ticks_diag + ticks_hsp + ticks_rsb;
	uint pct_seed = 0;
	uint pct_diag = 0;
	uint pct_hsp = 0;
	uint pct_rsb = 0;
	if (ticks_all > 0)
		{
		pct_seed = uint((ticks_seed * 100) / ticks_all);
		pct_diag = uint((ticks_diag * 100) / ticks_all);
		pct_hsp = uint((ticks_hsp * 100) / ticks_all);
		pct_rsb = uint((ticks_rsb * 100) / ticks_all);
		}
	ProgressLog("Kappa phases pct  seed=%u diag=%u hsp=%u rsb=%u\n",
		pct_seed, pct_diag, pct_hsp, pct_rsb);

	for (uint ti = 0; ti < nthreads; ++ti)
		{
		const KappaThreadDiag &d = m_thread_diags[ti];
		if (d.targets == 0 && d.batches == 0)
			continue;
		Log("  thr%u  batches=%u targets=%llu hsp=%llu hsp_cells=%llu secs=%u\n",
			ti, d.batches,
			(unsigned long long) d.targets,
			(unsigned long long) d.hsp,
			(unsigned long long) d.hsp_cells,
			d.secs);
		}
	}

void kappa_filter::run_filter(
	uint8_t **query_kappa_codeseqs,
	const uint *query_lengths,
	uint NQ,
	kappa_seqsource &db_ss)
	{
	m_query_kappa_codeseq_vec = query_kappa_codeseqs;
	m_query_lengths = query_lengths;
	m_QSeqCount = NQ;
	m_db_seqsource = &db_ss;
	reset_filter_stats();

	ProgressStep(0, 1000, "Kappa filter");
	time_t t_start = time(0);
	m_time_last_progress = t_start;
	RankedScoresBag::m_AnyLoScoreActive.store(false, std::memory_order_relaxed);

	const bool use_bcb_batch = (db_ss.m_KSSS == KSSS_bcb);
	vector<thread *> ts;
	uint ThreadCount = GetRequestedThreadCount();
	m_thread_diags.assign(ThreadCount, KappaThreadDiag());
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		if (use_bcb_batch)
			ts.push_back(new thread(static_bcb_thread_body, ThreadIndex));
		else
			ts.push_back(new thread(static_thread_body, ThreadIndex));
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	ProgressStep(999, 1000, "Kappa filter");

	uint total = kappa_filter::m_RSB.TruncateAllQueryVecs();

	uint rsb_sat = 0;
	const uint B = kappa_filter::m_RSB.m_B;
	for (uint q = 0; q < NQ; ++q)
		{
		if (SIZE(kappa_filter::m_RSB.m_QueryIdxToTargetIdxVec[q]) >= B)
			++rsb_sat;
		}

	uint NDB = 0;
	if (db_ss.m_bcb != 0)
		NDB = db_ss.m_bcb->GetChainCount();

	log_filter_diagnostics(NQ, NDB, total, rsb_sat);
	}
