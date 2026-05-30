#include "myutils.h"
#include "kappa_filter.h"
#include "kappa_mermx.h"
#include "kappa_dex.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "nu_filter.h"

void cmd_flat_search()
	{
	const string &QFN = g_Arg1;
	const string &DBFN = opt(db);

	BCAData QBCA;
	BCAData DBBCA;

	QBCA.Open(QFN);
	asserta(QBCA.m_HasNuSequences);

	flat_params params;
	params.init_from_cmdline();
	query_data **query_data_vec = QBCA.get_query_data_vec(params);

	DBBCA.Open(DBFN);

	const uint nquery = QBCA.GetChainCount();
	const uint TSeqCount = QBCA.GetChainCount();
	decide_query_or_db_kmer_neighborhood(nquery, TSeqCount);

	kappa_filter::init_kappa();
	kappa_filter::m_RSB.Init(nquery);

	kappa_dex QKmerIndex;
	QKmerIndex.Init();

	const uint k = flat_params::m_kappa_kmer_nrones;
	const kappa_mermx &GetKappaMerMx(uint k);
	const kappa_mermx &ScoreMx = GetKappaMerMx(k);
	asserta(ScoreMx.m_k == k);

	QKmerIndex.m_KmerSelfScores = ScoreMx.BuildSelfScores_Kmers();
	QKmerIndex.m_MinKmerSelfScore =  flat_params::m_kappa_min_mindiagscore;

	uint8_t **query_kappa_codeseqs = myalloc(uint8_t *, nquery);
	uint8_t **query_nu_codeseqs = myalloc(uint8_t *, nquery);
	uint *query_lengths = myalloc(uint, nquery);
	parasail_profile_t **query_parasail_profs = myalloc(parasail_profile_t *, nquery);
	parasail_profile_t **query_parasail_prof_revs = myalloc(parasail_profile_t *, nquery);

	for (uint chainidx = 0; chainidx < nquery; ++chainidx)
		{
		uint L = query_data_vec[chainidx]->m_chain->get_length();
		query_lengths[chainidx] = L;
		query_nu_codeseqs[chainidx] = query_data_vec[chainidx]->m_codeseq_nu;
		query_kappa_codeseqs[chainidx] = query_data_vec[chainidx]->m_codeseq_kappa;
		query_parasail_profs[chainidx] = query_data_vec[chainidx]->m_parasail_prof;
		query_parasail_prof_revs[chainidx] = query_data_vec[chainidx]->m_parasail_prof_rev;
		}

	QKmerIndex.from_codeseqs(
		query_kappa_codeseqs, query_lengths,
		QBCA.m_Labels, nquery);
#if DEBUG
	QKmerIndex.Validate();
#endif
	asserta(QKmerIndex.m_k == k);
	asserta(QKmerIndex.m_DictSize == flat_params::m_kappa_dict_size);
	asserta(ScoreMx.m_AS_pow[k] == QKmerIndex.m_DictSize);

	kappa_filter::m_ptrScoreMx = &ScoreMx;
	kappa_filter::m_ptrQKmerIndex = &QKmerIndex;

	kappa_seqsource db_ss;
	db_ss.OpenBCB(DBBCA);

	time_t t_kappa_filter_start = time(0);
	kappa_filter::run_filter(
		query_kappa_codeseqs, query_lengths, nquery, db_ss);

	vector<uint> dbidxs;
	unordered_map<uint, vector<uint> > dbidx_to_qidxs;
	kappa_filter::m_RSB.GetTargetInfo(dbidxs, dbidx_to_qidxs);
	time_t t_kappa_filter_end = time(0);
	uint kappa_filter_secs = uint(t_kappa_filter_end - t_kappa_filter_start);
	ProgressLog("Kappa filter %u secs\n", kappa_filter_secs);

	nu_filter::set_query_parasail_profiles(
		QBCA.m_Labels,
		query_parasail_profs,
		query_parasail_prof_revs,
		query_lengths, nquery);

	nu_filter::set_query_self_rev_scores(query_nu_codeseqs);
	nu_filter::run_filter(params, DBBCA, dbidxs, dbidx_to_qidxs);
	time_t t_nu_filter_end = time(0);
	uint nu_filter_secs = uint(t_nu_filter_end - t_kappa_filter_end);
	ProgressLog("Nu filter %u secs\n", kappa_filter_secs);
	}
