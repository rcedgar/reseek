#include "myutils.h"
#include "kappa_dex.h"
#include "kappa_mermx.h"
#include "kappa_filter.h"
#include "kappa_hsp.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "flat_chain.h"
#include "twohitdiag.h"
#include "rankedscoresbag.h"
#include "reseeker.h"
#include "bcadata.h"
#include <set>
#include <thread>
#include <atomic>

void set_default_stats();

// DB-index kappa prefilter (exact .kdx, hood on query) → HSP → RSB
// → optional search_post_kappa.
//   reseek -idx_search_kappa query.bcb -db db.kdx -input2 db.bcb
//     [-filter_only] [-output2 pairs.tsv] [-twohitdiag] [-output hits.tsv]
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

struct IdxSearchShared
	{
	const kappa_dex *Index = 0;
	const kappa_mermx *ScoreMx = 0;
	int MinScore = 0;
	uint8_t **q_kappa = 0;
	const uint *q_lengths = 0;
	const vector<string> *query_labels = 0;
	uint8_t **t_kappa = 0;
	const uint *t_lengths = 0;
	const vector<string> *t_labels = 0;
	uint nquery = 0;
	FILE *f_prehsp = 0;
	atomic<uint> next_qidx{0};
	atomic<uint64> n_prehsp{0};
	atomic<uint64> n_hsp_accept{0};
	atomic<time_t> time_last_progress{0};
	};

static void idx_search_thread_body(IdxSearchShared *S)
	{
	const kappa_dex &Index = *S->Index;
	const kappa_mermx &ScoreMx = *S->ScoreMx;
	const int MinScore = S->MinScore;
	const uint nseq = Index.m_nseq;

	uint *NeighborKmers = myalloc(uint, Index.m_DictSize);
	uint16_t *TBestScore = myalloc(uint16_t, nseq);
	uint32_t *THitList = myalloc(uint32_t, nseq);
	zero_array(TBestScore, nseq);

	TwoHitDiag Bag;
	set<uint32_t> OneHit;
	vector<uint> QKmers;
	vector<RankedScoreBatchEntry> Pending;
	Pending.reserve(kappa_filter::RSB_BATCH);

	uint64 local_prehsp = 0;
	uint64 local_hsp_accept = 0;
	uint counter = 0;

	for (;;)
		{
		const uint qidx = S->next_qidx.fetch_add(1, memory_order_relaxed);
		if (qidx >= S->nquery)
			break;

		if ((++counter) % 10 == 0)
			{
			time_t now = time(0);
			if (now > S->time_last_progress.load(memory_order_relaxed))
				{
				static mutex s_progress_lock;
				lock_guard<mutex> lock(s_progress_lock);
				const uint done = min(qidx + 1, S->nquery);
				uint pctx10 = (S->nquery == 0) ? 0 : (done * 1000u) / S->nquery;
				if (pctx10 >= 999) pctx10 = 998;
				ProgressStep(pctx10, 1000, "Kappa DB-index filter");
				S->time_last_progress.store(now, memory_order_relaxed);
				}
			}

		const uint QL = S->q_lengths[qidx];
		if (QL < flat_params::m_kappa_min_chainlength || QL < Index.m_K)
			continue;

		const byte *QSeq = S->q_kappa[qidx];
		const char *QLabel = (*S->query_labels)[qidx].c_str();

		Bag.Reset();
		OneHit.clear();
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
				const uint RowSize = Index.GetRowSize(Nbr);
				if (RowSize == 0)
					continue;
				uint DataOffset = Index.GetRowStart(Nbr);
				for (uint c = 0; c < RowSize; ++c)
					{
					uint32_t TSeqIdx;
					uint16_t TPos;
					Index.Get(DataOffset++, TSeqIdx, TPos);
					asserta(TSeqIdx < nseq);
					const uint16_t Diag = uint16_t(QL + TPos - QPos - 1);
					if (Diag > m_Mask14)
						continue;
					if (flat_params::m_kappa_onehitdiag)
						{
						asserta(TSeqIdx < UINT16_MAX);
						OneHit.insert((uint32_t(TSeqIdx) << 16) | uint32_t(Diag));
						}
					else
						Bag.Add(TSeqIdx, Diag);
					}
				}
			}

		uint n_thit = 0;

		auto NoteBest = [&](uint tidx, uint16_t diag)
			{
			++local_prehsp;
			if (S->f_prehsp != 0)
				{
				lock_guard<mutex> lock(kappa_filter::m_prehsp_dump_mutex);
				kappa_filter::write_prehsp_hit(S->f_prehsp, QLabel,
					(*S->t_labels)[tidx].c_str(), qidx, tidx, diag);
				}

			const uint TL = S->t_lengths[tidx];
			if (TL < flat_params::m_kappa_min_chainlength)
				return;

			int DiagScore = ExtendDiagToHSP_DB(QSeq, QL, S->t_kappa[tidx], TL,
				diag, qidx, kappa_filter::m_RSB);
			if (DiagScore <= 0)
				return;
			if (DiagScore < flat_params::m_kappa_min_diagscore)
				return;
			if (DiagScore >= UINT16_MAX)
				DiagScore = UINT16_MAX - 1;
			const uint16_t sc = uint16_t(DiagScore);
			if (TBestScore[tidx] == 0)
				THitList[n_thit++] = tidx;
			if (sc > TBestScore[tidx])
				TBestScore[tidx] = sc;
			};

		if (flat_params::m_kappa_onehitdiag)
			{
			for (set<uint32_t>::const_iterator iter = OneHit.begin();
				 iter != OneHit.end(); ++iter)
				{
				const uint32_t pair = *iter;
				NoteBest(pair >> 16, uint16_t(pair & 0xffff));
				}
			}
		else
			{
			if (flat_params::m_kappa_twohitdiag)
				Bag.SetDupes();
			else
				Bag.SetUniqueFine();
			for (uint i = 0; i < Bag.m_DupeCount; ++i)
				NoteBest(Bag.m_DupeSeqIdxs[i], Bag.m_DupeDiags[i]);
			}

		for (uint i = 0; i < n_thit; ++i)
			{
			const uint tidx = THitList[i];
			const uint16_t sc = TBestScore[tidx];
			TBestScore[tidx] = 0;
			RankedScoreBatchEntry e;
			e.QueryIdx = qidx;
			e.TargetIdx = tidx;
			e.Score = sc;
			Pending.push_back(e);
			++local_hsp_accept;
			if (Pending.size() >= kappa_filter::RSB_BATCH)
				kappa_filter::m_RSB.AddScoresBatch(Pending);
			}
		}

	if (!Pending.empty())
		kappa_filter::m_RSB.AddScoresBatch(Pending);

	S->n_prehsp.fetch_add(local_prehsp, memory_order_relaxed);
	S->n_hsp_accept.fetch_add(local_hsp_accept, memory_order_relaxed);

	myfree(NeighborKmers);
	myfree(TBestScore);
	myfree(THitList);
	}

void cmd_idx_search_kappa()
	{
	asserta(optset_db);
	asserta(optset_input2);
	asserta(EndsWith(g_Arg1, ".bcb"));
	asserta(EndsWith(string(opt(input2)), ".bcb"));

	set_default_stats();
	flat_params params;
	params.init_from_cmdline();
	params.logme();

	g_flat_n_truncated_chains = 0;

	kappa_filter::init_kappa();

	kappa_dex Index;
	Index.FromFile(opt(db));
	asserta(Index.m_nseq > 0);
	asserta(Index.m_DictSize > 0);

	const uint k = Index.m_k;
	const kappa_mermx &GetKappaMerMx(uint k);
	const kappa_mermx &ScoreMx = GetKappaMerMx(k);
	asserta(ScoreMx.m_k == k);
	asserta(ScoreMx.m_AS_pow[k] == Index.m_DictSize);

	Index.m_KmerSelfScores = ScoreMx.BuildSelfScores_Kmers();
	const int MinScore = flat_params::m_kappa_min_kmerpairscore;
	if (MinScore != Index.m_MinKmerSelfScore)
		ProgressLog("Warning: -kappa_minkmerscore %d != index MinKmerSelfScore %d\n",
			MinScore, Index.m_MinKmerSelfScore);

	if (flat_params::m_kappa_onehitdiag && Index.m_nseq >= UINT16_MAX)
		Die("-onehitdiag idx_search_kappa requires db nseq < 65535");

	ProgressLog("Index db k-mer neighborhoods (exact .kdx, hood on query)\n");

	BCAData DBBCA;
	DBBCA.Open(opt(input2));
	asserta(DBBCA.m_HasNuSequences);
	asserta(DBBCA.GetChainCount() == Index.m_nseq);

	uint8_t **t_kappa = 0;
	uint *t_lengths = 0;
	DBBCA.make_kappa_codeseqs(&t_kappa, &t_lengths);
	const vector<string> &t_labels = DBBCA.m_Labels;

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
		fprintf(f_prehsp, "# index\t%s\n", opt(db));
		fprintf(f_prehsp, "# query\t%s\n", g_Arg1.c_str());
		fprintf(f_prehsp, "# targets\t%s\n", opt(input2));
		}

	IdxSearchShared Shared;
	Shared.Index = &Index;
	Shared.ScoreMx = &ScoreMx;
	Shared.MinScore = MinScore;
	Shared.q_kappa = q_kappa;
	Shared.q_lengths = q_lengths;
	Shared.query_labels = &query_labels;
	Shared.t_kappa = t_kappa;
	Shared.t_lengths = t_lengths;
	Shared.t_labels = &t_labels;
	Shared.nquery = nquery;
	Shared.f_prehsp = f_prehsp;
	Shared.time_last_progress = time(0);

	const uint ThreadCount = GetRequestedThreadCount();
	ProgressLog("Kappa DB-index filter threads %u  queries %u\n",
		ThreadCount, nquery);
	ProgressStep(0, 1000, "Kappa DB-index filter");
	time_t t0 = time(0);

	vector<thread *> ts;
	for (uint ti = 0; ti < ThreadCount; ++ti)
		ts.push_back(new thread(idx_search_thread_body, &Shared));
	for (uint ti = 0; ti < ThreadCount; ++ti)
		ts[ti]->join();
	for (uint ti = 0; ti < ThreadCount; ++ti)
		delete ts[ti];

	ProgressStep(999, 1000, "Kappa DB-index filter");

	if (f_prehsp != 0)
		CloseStdioFile(f_prehsp);

	uint total = kappa_filter::m_RSB.TruncateAllQueryVecs();
	time_t t1 = time(0);
	ProgressLog("Kappa DB-index filter %u secs  prehsp=%llu  rsb_pairs=%u  hsp_targets=%llu\n",
		uint(t1 - t0),
		(unsigned long long) Shared.n_prehsp.load(),
		total,
		(unsigned long long) Shared.n_hsp_accept.load());

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
	}
