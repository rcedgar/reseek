#include "myutils.h"
#include "kappa_dex.h"
#include "kappa_mermx.h"
#include "kappa_filter.h"
#include "kappa_hsp.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "flat_chain.h"
#include "rankedscoresbag.h"
#include "reseeker.h"
#include "bcadata.h"
#include <algorithm>
#include <thread>
#include <atomic>
#include <unordered_set>

void set_default_stats();

static const uint32_t IDX_DIAG_MASK14 = 0b11111111111111;
// Fixed per-thread raw-seed cap (Foldseek-style). Do NOT scale with nseq —
// 2*nseq on ~50M DBs yields ~1e8 and stalls inside a single query.
static const uint IDX_SEED_BUF_CAP = 1000000;

static int ExtendDiagToHSP_DB(
	const byte *QSeq, uint QL,
	const byte *TSeq, uint TL,
	uint16_t Diag,
	uint QSeqIdx,
	RankedScoresBag &RSB)
	{
	if (flat_params::m_kappa_hsp_rsb_prune &&
		RSB.m_AnyLoScoreActive.load(std::memory_order_relaxed))
		{
		const uint16_t LoScore = RSB.GetLoScore(QSeqIdx);
		if (LoScore > 0)
			{
			int mini, minj, n;
			kappa_get_hsp_limits(int(QL), int(TL), int(Diag),
				mini, minj, n);
			const int bound = n * flat_params::m_kappa_max_pos_logodds;
			if (bound < LoScore)
				return 0;
			}
		}
	return kappa_find_hsp(QSeq, TSeq, int(QL), int(TL), int(Diag));
	}

static inline uint64_t PackTDiag(uint32_t tidx, uint16_t diag)
	{
	return (uint64_t(tidx) << 16) | uint64_t(diag);
	}

struct IdxSeed
	{
	uint32_t tidx;
	uint16_t diag;
	};

struct IdxSearchShared
	{
	const kappa_mermx *ScoreMx = 0;
	int16_t *KmerSelfScores = 0;	// shared; used by shard builds
	int MinScore = 0;
	uint8_t **q_kappa = 0;
	const uint *q_lengths = 0;
	const vector<string> *query_labels = 0;
	uint8_t **t_kappa = 0;
	const uint *t_lengths = 0;
	const vector<string> *t_labels = 0;
	uint nquery = 0;
	uint progress_total = 0;	// units for ProgressStep (nquery, or shards*nquery)
	uint seed_cap = 0;
	uint max_seqs = 0;	// 0 = unlimited; else Foldseek-like top-N HSP floor
	FILE *f_prehsp = 0;
	atomic<uint> next_qidx{0};
	atomic<uint> n_queries_done{0};
	atomic<uint64> n_prehsp{0};
	atomic<uint64> n_hsp_accept{0};
	atomic<uint64> n_hsp_skip{0};
	atomic<uint64> n_seed_flushes{0};
	atomic<time_t> time_last_progress{0};
	};

struct IdxThreadCtx
	{
	IdxSearchShared *S = 0;
	const kappa_dex *Index = 0;
	uint seq_base = 0;	// global target index = seq_base + local index id
	bool all_queries = false;	// true: stream every query (shard mode)
	};

static void MaybeReportIdxProgress(IdxSearchShared *S)
	{
	const time_t now = time(0);
	if (now <= S->time_last_progress.load(memory_order_relaxed))
		return;
	static mutex s_progress_lock;
	lock_guard<mutex> lock(s_progress_lock);
	if (now <= S->time_last_progress.load(memory_order_relaxed))
		return;
	const uint ntotal = (S->progress_total > 0) ? S->progress_total : S->nquery;
	if (ntotal == 0)
		return;
	uint done = S->n_queries_done.load(memory_order_relaxed);
	if (done > ntotal)
		done = ntotal;
	// ProgressStep(i,N) requires i < N; i==0 resets internals — avoid after start.
	uint i = done;
	if (i >= ntotal)
		i = ntotal - 1;
	if (i == 0)
		i = 1;
	ProgressStep(i, ntotal,
		"Kappa DB-index filter  done=%u/%u  flushes=%llu  prehsp=%llu",
		done, ntotal,
		(unsigned long long) S->n_seed_flushes.load(memory_order_relaxed),
		(unsigned long long) S->n_prehsp.load(memory_order_relaxed));
	S->time_last_progress.store(now, memory_order_relaxed);
	}

static void idx_search_thread_body(IdxThreadCtx *Ctx)
	{
	IdxSearchShared *S = Ctx->S;
	const kappa_dex &Index = *Ctx->Index;
	const kappa_mermx &ScoreMx = *S->ScoreMx;
	const int MinScore = S->MinScore;
	const uint nseq = Index.m_nseq;
	const uint SeedCap = S->seed_cap;
	const uint MaxSeqs = S->max_seqs;
	const bool twohit = flat_params::m_kappa_twohitdiag;
	const uint seq_base = Ctx->seq_base;
	const bool all_queries = Ctx->all_queries;

	uint *NeighborKmers = myalloc(uint, Index.m_DictSize);
	uint16_t *TBestScore = myalloc(uint16_t, nseq);
	uint32_t *THitList = myalloc(uint32_t, nseq);
	zero_array(TBestScore, nseq);

	vector<IdxSeed> SeedBuf;
	SeedBuf.reserve(SeedCap);
	unordered_set<uint64_t> TwoHitCarry;
	if (twohit)
		TwoHitCarry.reserve(SeedCap);

	vector<uint> QKmers;
	vector<RankedScoreBatchEntry> Pending;
	Pending.reserve(kappa_filter::RSB_BATCH);

	uint64 local_prehsp = 0;
	uint64 local_hsp_accept = 0;
	uint64 local_hsp_skip = 0;
	uint n_thit = 0;
	uint16_t max_seqs_floor = 0;

	auto RecomputeMaxSeqsFloor = [&]()
		{
		if (MaxSeqs == 0 || n_thit < MaxSeqs)
			{
			max_seqs_floor = 0;
			return;
			}
		uint16_t mn = UINT16_MAX;
		for (uint i = 0; i < n_thit; ++i)
			{
			const uint16_t sc = TBestScore[THitList[i]];
			if (sc < mn)
				mn = sc;
			}
		max_seqs_floor = (mn == UINT16_MAX) ? 0 : mn;
		};

	auto EvictWorstIfNeeded = [&]()
		{
		if (MaxSeqs == 0)
			return;
		while (n_thit > MaxSeqs)
			{
			uint worst_i = 0;
			uint16_t worst_sc = TBestScore[THitList[0]];
			for (uint i = 1; i < n_thit; ++i)
				{
				const uint16_t sc = TBestScore[THitList[i]];
				if (sc < worst_sc)
					{
					worst_sc = sc;
					worst_i = i;
					}
				}
			TBestScore[THitList[worst_i]] = 0;
			THitList[worst_i] = THitList[--n_thit];
			}
		RecomputeMaxSeqsFloor();
		};

	auto NoteBest = [&](uint qidx, const byte *QSeq, uint QL,
		const char *QLabel, uint tidx_local, uint16_t diag)
		{
		const uint tidx_global = seq_base + tidx_local;
		++local_prehsp;
		if (S->f_prehsp != 0)
			{
			lock_guard<mutex> lock(kappa_filter::m_prehsp_dump_mutex);
			kappa_filter::write_prehsp_hit(S->f_prehsp, QLabel,
				(*S->t_labels)[tidx_global].c_str(), qidx, tidx_global, diag);
			}

		const uint TL = S->t_lengths[tidx_global];
		if (TL < flat_params::m_kappa_min_chainlength)
			return;

		const bool already = (TBestScore[tidx_local] != 0);

		// Foldseek-like: once top-N is full, skip HSP that cannot beat the floor.
		if (!already && max_seqs_floor > 0)
			{
			int mini, minj, n;
			kappa_get_hsp_limits(int(QL), int(TL), int(diag), mini, minj, n);
			const int bound = n * flat_params::m_kappa_max_pos_logodds;
			if (bound < int(max_seqs_floor))
				{
				++local_hsp_skip;
				return;
				}
			}

		int DiagScore = ExtendDiagToHSP_DB(QSeq, QL, S->t_kappa[tidx_global], TL,
			diag, qidx, kappa_filter::m_RSB);
		if (DiagScore <= 0)
			return;
		if (DiagScore < flat_params::m_kappa_min_diagscore)
			return;
		if (DiagScore >= UINT16_MAX)
			DiagScore = UINT16_MAX - 1;
		const uint16_t sc = uint16_t(DiagScore);

		if (!already && max_seqs_floor > 0 && sc <= max_seqs_floor)
			{
			++local_hsp_skip;
			return;
			}

		if (!already)
			THitList[n_thit++] = tidx_local;
		if (sc > TBestScore[tidx_local])
			TBestScore[tidx_local] = sc;
		EvictWorstIfNeeded();
		};

	auto FlushSeeds = [&](uint qidx, const byte *QSeq, uint QL, const char *QLabel)
		{
		const uint n = SIZE(SeedBuf);
		if (n == 0)
			return;
		S->n_seed_flushes.fetch_add(1, memory_order_relaxed);

		sort(SeedBuf.begin(), SeedBuf.end(),
			[](const IdxSeed &a, const IdxSeed &b)
				{
				if (a.tidx != b.tidx)
					return a.tidx < b.tidx;
				return a.diag < b.diag;
				});

		for (uint i = 0; i < n; )
			{
			const uint32_t tidx = SeedBuf[i].tidx;
			const uint16_t diag = SeedBuf[i].diag;
			uint j = i + 1;
			while (j < n && SeedBuf[j].tidx == tidx && SeedBuf[j].diag == diag)
				++j;
			const uint cnt = j - i;
			i = j;

			if (twohit)
				{
				if (cnt >= 2)
					{
					TwoHitCarry.erase(PackTDiag(tidx, diag));
					NoteBest(qidx, QSeq, QL, QLabel, tidx, diag);
					}
				else
					{
					const uint64_t key = PackTDiag(tidx, diag);
					if (TwoHitCarry.find(key) != TwoHitCarry.end())
						{
						TwoHitCarry.erase(key);
						NoteBest(qidx, QSeq, QL, QLabel, tidx, diag);
						}
					else if (TwoHitCarry.size() < SeedCap)
						TwoHitCarry.insert(key);
					}
				}
			else
				NoteBest(qidx, QSeq, QL, QLabel, tidx, diag);
			}

		SeedBuf.clear();
		if (local_prehsp > 0)
			{
			S->n_prehsp.fetch_add(local_prehsp, memory_order_relaxed);
			local_prehsp = 0;
			}
		MaybeReportIdxProgress(S);
		};

	auto PushSeed = [&](uint qidx, const byte *QSeq, uint QL, const char *QLabel,
		uint32_t tidx, uint16_t diag)
		{
		if (SIZE(SeedBuf) >= SeedCap)
			FlushSeeds(qidx, QSeq, QL, QLabel);
		IdxSeed s;
		s.tidx = tidx;
		s.diag = diag;
		SeedBuf.push_back(s);
		};

	auto ProcessOneQuery = [&](uint qidx)
		{
		MaybeReportIdxProgress(S);

		const uint QL = S->q_lengths[qidx];
		if (QL < flat_params::m_kappa_min_chainlength || QL < Index.m_K)
			{
			S->n_queries_done.fetch_add(1, memory_order_relaxed);
			return;
			}

		const byte *QSeq = S->q_kappa[qidx];
		const char *QLabel = (*S->query_labels)[qidx].c_str();

		SeedBuf.clear();
		TwoHitCarry.clear();
		n_thit = 0;
		max_seqs_floor = 0;

		Index.GetKmers(QSeq, QL, QKmers);
		const uint NK = SIZE(QKmers);

		for (uint QPos = 0; QPos < NK; ++QPos)
			{
			uint QKmer = QKmers[QPos];
			if (QKmer == UINT_MAX)
				continue;
			asserta(QKmer < Index.m_DictSize);
			if (Index.m_KmerSelfScores[QKmer] < MinScore)
				continue;

			const uint n = ScoreMx.GetHighScoringKmers(QKmer, short(MinScore),
				NeighborKmers);
			for (uint j = 0; j < n; ++j)
				{
				const uint Nbr = NeighborKmers[j];
				const uint64_t RowSize = Index.GetRowSize(Nbr);
				if (RowSize == 0)
					continue;
				uint64_t DataOffset = Index.GetRowStart(Nbr);
				for (uint64_t c = 0; c < RowSize; ++c)
					{
					uint32_t TSeqIdx;
					uint16_t TPos;
					Index.Get(DataOffset++, TSeqIdx, TPos);
					asserta(TSeqIdx < nseq);
					const uint16_t Diag = uint16_t(QL + TPos - QPos - 1);
					if (Diag > IDX_DIAG_MASK14)
						continue;
					PushSeed(qidx, QSeq, QL, QLabel, TSeqIdx, Diag);
					}
				}
			}

		FlushSeeds(qidx, QSeq, QL, QLabel);
		TwoHitCarry.clear();

		for (uint i = 0; i < n_thit; ++i)
			{
			const uint tidx_local = THitList[i];
			const uint16_t sc = TBestScore[tidx_local];
			TBestScore[tidx_local] = 0;
			RankedScoreBatchEntry e;
			e.QueryIdx = qidx;
			e.TargetIdx = seq_base + tidx_local;
			e.Score = sc;
			Pending.push_back(e);
			++local_hsp_accept;
			if (Pending.size() >= kappa_filter::RSB_BATCH)
				kappa_filter::m_RSB.AddScoresBatch(Pending);
			}

		S->n_queries_done.fetch_add(1, memory_order_relaxed);
		if (local_prehsp > 0)
			{
			S->n_prehsp.fetch_add(local_prehsp, memory_order_relaxed);
			local_prehsp = 0;
			}
		MaybeReportIdxProgress(S);
		};

	if (all_queries)
		{
		for (uint qidx = 0; qidx < S->nquery; ++qidx)
			ProcessOneQuery(qidx);
		}
	else
		{
		for (;;)
			{
			const uint qidx = S->next_qidx.fetch_add(1, memory_order_relaxed);
			if (qidx >= S->nquery)
				break;
			ProcessOneQuery(qidx);
			}
		}

	if (!Pending.empty())
		kappa_filter::m_RSB.AddScoresBatch(Pending);

	S->n_prehsp.fetch_add(local_prehsp, memory_order_relaxed);
	S->n_hsp_accept.fetch_add(local_hsp_accept, memory_order_relaxed);
	S->n_hsp_skip.fetch_add(local_hsp_skip, memory_order_relaxed);

	myfree(NeighborKmers);
	myfree(TBestScore);
	myfree(THitList);
	}

struct IdxShardArgs
	{
	IdxSearchShared *S = 0;
	uint seq_lo = 0;
	uint seq_hi = 0;
	uint shard_id = 0;
	uint nshards = 0;
	};

static void idx_search_shard_thread(IdxShardArgs *A)
	{
	IdxSearchShared *S = A->S;
	const uint lo = A->seq_lo;
	const uint hi = A->seq_hi;
	asserta(hi >= lo);
	const uint nshard = hi - lo;
	if (nshard == 0)
		{
		// Still count queries as done so progress reaches 100%.
		S->n_queries_done.fetch_add(S->nquery, memory_order_relaxed);
		return;
		}

	kappa_dex Index;
	Index.Init();
	asserta(S->KmerSelfScores != 0);
	Index.m_KmerSelfScores = S->KmerSelfScores;
	Index.m_MinKmerSelfScore = S->MinScore;
	Index.m_AddNeighborhood = false;
	Index.m_UniqueKmer = opt(unique_kmer);
	Index.m_ptrScoreMx = 0;
	// Serial build: shard workers already run in parallel.
	Index.from_codeseqs(S->t_kappa + lo, S->t_lengths + lo, *S->t_labels,
		nshard, /*build_threads=*/1);

	IdxThreadCtx Ctx;
	Ctx.S = S;
	Ctx.Index = &Index;
	Ctx.seq_base = lo;
	Ctx.all_queries = true;
	idx_search_thread_body(&Ctx);
	}

void cmd_idx_search_kappa()
	{
	asserta(optset_db);
	asserta(EndsWith(g_Arg1, ".bcb"));
	asserta(EndsWith(string(opt(db)), ".bcb"));

	set_default_stats();
	flat_params params;
	params.init_from_cmdline();
	params.logme();

	g_flat_n_truncated_chains = 0;

	kappa_filter::init_kappa();

	const uint DbShards = optset_db_shards ? opt(db_shards) : 1;
	if (DbShards == 0)
		Die("-db_shards must be >= 1 (omit or use 1 for unsharded path)");
	if (DbShards > 1 && optset_kdx)
		Die("-db_shards > 1 builds per-shard indexes; do not pass -kdx");

	BCAData DBBCA;
	DBBCA.Open(opt(db));
	asserta(DBBCA.m_HasNuSequences);

	uint8_t **t_kappa = 0;
	uint *t_lengths = 0;
	DBBCA.make_kappa_codeseqs(&t_kappa, &t_lengths);
	const vector<string> &t_labels = DBBCA.m_Labels;
	const uint nseq = DBBCA.GetChainCount();

	kappa_dex Index;
	const kappa_mermx &GetKappaMerMx(uint k);
	const int MinScore = flat_params::m_kappa_min_kmerpairscore;
	const kappa_mermx *ptrScoreMx = 0;

	if (DbShards > 1)
		{
		const uint k = flat_params::m_kappa_kmer_nrones;
		ptrScoreMx = &GetKappaMerMx(k);
		asserta(ptrScoreMx->m_k == k);
		ProgressLog("DB-shard mode  shards=%u  nseq=%u  (private index per shard, all queries)\n",
			DbShards, nseq);
		}
	else if (optset_kdx)
		{
		Index.FromFile(opt(kdx));
		asserta(Index.m_nseq > 0);
		asserta(Index.m_DictSize > 0);
		asserta(Index.m_nseq == nseq);
		const uint k = Index.m_k;
		ptrScoreMx = &GetKappaMerMx(k);
		asserta(ptrScoreMx->m_k == k);
		asserta(ptrScoreMx->m_AS_pow[k] == Index.m_DictSize);
		Index.m_KmerSelfScores = ptrScoreMx->BuildSelfScores_Kmers();
		if (MinScore != Index.m_MinKmerSelfScore)
			ProgressLog("Warning: -kappa_minkmerscore %d != index MinKmerSelfScore %d\n",
				MinScore, Index.m_MinKmerSelfScore);
		ProgressLog("Loaded .kdx %s  (exact DB index, hood on query)\n", opt(kdx));
		}
	else
		{
		Index.Init();
		const uint k = flat_params::m_kappa_kmer_nrones;
		ptrScoreMx = &GetKappaMerMx(k);
		asserta(ptrScoreMx->m_k == k);
		Index.m_KmerSelfScores = ptrScoreMx->BuildSelfScores_Kmers();
		Index.m_MinKmerSelfScore = MinScore;
		Index.m_AddNeighborhood = false;
		Index.m_UniqueKmer = opt(unique_kmer);
		Index.m_ptrScoreMx = 0;
		ProgressLog("Building kappa_dex on the fly from %s%s\n",
			opt(db), Index.m_UniqueKmer ? "  (-unique_kmer)" : "");
		Index.from_codeseqs(t_kappa, t_lengths, t_labels, nseq);
		asserta(Index.m_nseq == nseq);
		ProgressLog("On-the-fly index ready  nseq=%u  postings=%s\n",
			nseq, Int64ToStr(Index.m_Size));
		}

	const kappa_mermx &ScoreMx = *ptrScoreMx;

	BCAData QBCA;
	QBCA.Open(g_Arg1);
	asserta(QBCA.m_HasNuSequences);

	vector<string> query_labels;
	struct_data **struct_data_vec = 0;
	uint8_t **q_kappa = 0;
	uint *q_lengths = 0;
	uint nquery = 0;

	const flat_chain_t **query_chains = 0;
	uint8_t **query_nu_codeseqs = 0;
	uint8_t **query_mega_profs = 0;
	const float **query_mega_pssms = 0;
	const float **query_mega_pssm_revs = 0;
	const sid_t **query_distmxs = 0;
	parasail_profile_t **query_parasail_profs = 0;
	parasail_profile_t **query_parasail_prof_revs = 0;

	if (opt(filter_only))
		{
		QBCA.make_kappa_codeseqs(&q_kappa, &q_lengths);
		nquery = QBCA.GetChainCount();
		query_labels = QBCA.m_Labels;
		}
	else
		{
		struct_data_vec = QBCA.get_struct_data_vec(params, query_labels);
		nquery = uint(query_labels.size());
		q_kappa = myalloc(uint8_t *, nquery);
		q_lengths = myalloc(uint, nquery);
		query_chains = myalloc(const flat_chain_t *, nquery);
		query_nu_codeseqs = myalloc(uint8_t *, nquery);
		query_mega_profs = myalloc(uint8_t *, nquery);
		query_mega_pssms = myalloc(const float *, nquery);
		query_mega_pssm_revs = myalloc(const float *, nquery);
		query_distmxs = myalloc(const sid_t *, nquery);
		query_parasail_profs = myalloc(parasail_profile_t *, nquery);
		query_parasail_prof_revs = myalloc(parasail_profile_t *, nquery);
		for (uint i = 0; i < nquery; ++i)
			{
			q_lengths[i] = struct_data_vec[i]->m_chain->get_length();
			q_kappa[i] = struct_data_vec[i]->m_codeseq_kappa;
			query_chains[i] = struct_data_vec[i]->m_chain;
			query_nu_codeseqs[i] = struct_data_vec[i]->m_codeseq_nu;
			query_mega_profs[i] = struct_data_vec[i]->m_mega_prof;
			query_mega_pssms[i] = struct_data_vec[i]->m_mega_pssm;
			query_mega_pssm_revs[i] = struct_data_vec[i]->m_mega_pssm_rev;
			query_distmxs[i] = struct_data_vec[i]->m_distmx;
			query_parasail_profs[i] = struct_data_vec[i]->m_parasail_prof;
			query_parasail_prof_revs[i] = struct_data_vec[i]->m_parasail_prof_rev;
			}
		}

	kappa_filter::m_RSB.Init(nquery);
	RankedScoresBag::m_AnyLoScoreActive.store(false, std::memory_order_relaxed);

	FILE *f_prehsp = 0;
	if (optset_dump_prefilter_prehsp)
		{
		f_prehsp = CreateStdioFile(opt(dump_prefilter_prehsp));
		kappa_filter::write_prehsp_tsv_header(f_prehsp, "db_index");
		if (DbShards > 1)
			fprintf(f_prehsp, "# index\tdb_shards\t%u\n", DbShards);
		else if (optset_kdx)
			fprintf(f_prehsp, "# index\t%s\n", opt(kdx));
		else
			fprintf(f_prehsp, "# index\ton_the_fly\n");
		fprintf(f_prehsp, "# query\t%s\n", g_Arg1.c_str());
		fprintf(f_prehsp, "# targets\t%s\n", opt(db));
		}

	const uint SeedCap = IDX_SEED_BUF_CAP;

	IdxSearchShared Shared;
	Shared.ScoreMx = &ScoreMx;
	Shared.MinScore = MinScore;
	Shared.q_kappa = q_kappa;
	Shared.q_lengths = q_lengths;
	Shared.query_labels = &query_labels;
	Shared.t_kappa = t_kappa;
	Shared.t_lengths = t_lengths;
	Shared.t_labels = &t_labels;
	Shared.nquery = nquery;
	Shared.seed_cap = SeedCap;
	Shared.max_seqs = optset_max_seqs ? opt(max_seqs) : 0;
	Shared.f_prehsp = f_prehsp;
	Shared.time_last_progress = time(0);

	int16_t *ShardSelfScores = 0;
	if (DbShards > 1)
		{
		ShardSelfScores = ScoreMx.BuildSelfScores_Kmers();
		Shared.KmerSelfScores = ShardSelfScores;
		}

	asserta(nquery > 0);
	time_t t0 = time(0);

	if (DbShards > 1)
		{
		uint T = DbShards;
		if (T > nseq)
			T = nseq;
		Shared.progress_total = T * nquery;
		if (Shared.max_seqs > 0)
			ProgressLog("Kappa DB-index filter db_shards %u  queries %u  seed_buf %u  max_seqs %u\n",
				T, nquery, SeedCap, Shared.max_seqs);
		else
			ProgressLog("Kappa DB-index filter db_shards %u  queries %u  seed_buf %u  max_seqs=off\n",
				T, nquery, SeedCap);
		ProgressStep(0, Shared.progress_total, "Kappa DB-index filter");

		vector<IdxShardArgs> args(T);
		vector<thread *> ts;
		for (uint ti = 0; ti < T; ++ti)
			{
			const uint lo = uint((uint64_t(ti) * nseq) / T);
			const uint hi = uint((uint64_t(ti + 1) * nseq) / T);
			args[ti].S = &Shared;
			args[ti].seq_lo = lo;
			args[ti].seq_hi = hi;
			args[ti].shard_id = ti;
			args[ti].nshards = T;
			ts.push_back(new thread(idx_search_shard_thread, &args[ti]));
			}
		for (uint ti = 0; ti < T; ++ti)
			ts[ti]->join();
		for (uint ti = 0; ti < T; ++ti)
			delete ts[ti];

		ProgressStep(Shared.progress_total - 1, Shared.progress_total,
			"Kappa DB-index filter");
		}
	else
		{
		const uint ThreadCount = GetRequestedThreadCount();
		Shared.progress_total = nquery;
		if (Shared.max_seqs > 0)
			ProgressLog("Kappa DB-index filter threads %u  queries %u  seed_buf %u  max_seqs %u\n",
				ThreadCount, nquery, SeedCap, Shared.max_seqs);
		else
			ProgressLog("Kappa DB-index filter threads %u  queries %u  seed_buf %u  max_seqs=off\n",
				ThreadCount, nquery, SeedCap);
		ProgressStep(0, nquery, "Kappa DB-index filter");

		IdxThreadCtx Ctx;
		Ctx.S = &Shared;
		Ctx.Index = &Index;
		Ctx.seq_base = 0;
		Ctx.all_queries = false;

		vector<thread *> ts;
		for (uint ti = 0; ti < ThreadCount; ++ti)
			ts.push_back(new thread(idx_search_thread_body, &Ctx));
		for (uint ti = 0; ti < ThreadCount; ++ti)
			ts[ti]->join();
		for (uint ti = 0; ti < ThreadCount; ++ti)
			delete ts[ti];

		ProgressStep(nquery - 1, nquery, "Kappa DB-index filter");
		}

	if (f_prehsp != 0)
		CloseStdioFile(f_prehsp);

	uint total = kappa_filter::m_RSB.TruncateAllQueryVecs();
	time_t t1 = time(0);
	ProgressLog("Kappa DB-index filter %u secs  prehsp=%s  rsb_pairs=%u  hsp_targets=%llu  hsp_skip=%llu  flushes=%llu\n",
		uint(t1 - t0),
		FloatToStr((float) Shared.n_prehsp.load()),
		total,
		(unsigned long long) Shared.n_hsp_accept.load(),
		(unsigned long long) Shared.n_hsp_skip.load(),
		(unsigned long long) Shared.n_seed_flushes.load());

	vector<uint> dbidxs;
	unordered_map<uint, vector<uint> > dbidx_to_qidxs;
	unordered_map<uint, vector<uint> > dbidx_to_diagscores;
	uint max_queries_per_target = 0;
	kappa_filter::m_RSB.GetTargetInfoSorted(
		dbidxs, dbidx_to_qidxs, dbidx_to_diagscores,
		max_queries_per_target);

	if (optset_output2)
		{
		const string &fn = opt(output2);
		FILE *f = CreateStdioFile(fn);
		Progress("Writing %s ...", fn.c_str());
		for (size_t i = 0; i < dbidxs.size(); ++i)
			{
			uint dbidx = dbidxs[i];
			const vector<uint> &qidxs = dbidx_to_qidxs[dbidx];
			for (auto qidx : qidxs)
				{
				fprintf(f, "%s\t%s\n",
					query_labels[qidx].c_str(),
					t_labels[dbidx].c_str());
				}
			}
		Progress(" done\n");
		CloseStdioFile(f);
		}

	if (opt(filter_only))
		{
		QBCA.Close();
		DBBCA.Close();
		return;
		}

	reseeker::set_params(params);
	reseeker::set_query_data(
		query_chains,
		query_labels,
		query_parasail_profs,
		query_parasail_prof_revs,
		query_mega_pssms,
		query_mega_pssm_revs,
		query_distmxs,
		q_lengths,
		nquery);

	if (optset_mints) reseeker::m_mints = float(opt(mints));
	reseeker::set_query_self_rev_scores(query_nu_codeseqs);
	reseeker::set_query_mega_self_rev_scores(query_mega_profs);
	reseeker::m_max_queries_per_target = max_queries_per_target;
	reseeker::search_post_kappa(DBBCA, dbidxs, dbidx_to_qidxs, dbidx_to_diagscores);
	time_t t2 = time(0);
	ProgressLog("Nu filter %u secs\n", uint(t2 - t1));

	QBCA.Close();
	DBBCA.Close();
	if (ShardSelfScores != 0)
		myfree(ShardSelfScores);
	}
