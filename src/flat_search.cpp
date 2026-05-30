#include "myutils.h"
#include "kappa_filter.h"
#include "kappa_mermx.h"
#include "kappa_dex.h"
#include "flat_params.h"
#include "flat_helpers.h"

void cmd_flat_search()
	{
	const string &QFN = g_Arg1;
	const string &DBFN = opt(db);

	BCAData QBCA;
	BCAData TBCA;

	QBCA.Open(QFN);
	asserta(QBCA.m_HasNuSequences);

	flat_params params;
	params.init_from_cmdline();
	query_data **query_data_vec = QBCA.get_query_data_vec(params);

	TBCA.Open(DBFN);

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

	uint8_t **query_kappa_codeseqs2 = myalloc(uint8_t *, nquery);
	uint *query_lengths2 = myalloc(uint, nquery);
	for (uint chainidx = 0; chainidx < nquery; ++chainidx)
		{
		uint L = query_data_vec[chainidx]->m_chain->get_length();
		query_lengths2[chainidx] = L;
		query_kappa_codeseqs2[chainidx] = query_data_vec[chainidx]->m_codeseq_kappa;
		}

	QKmerIndex.from_codeseqs(
		query_kappa_codeseqs2, query_lengths2,
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
	db_ss.OpenBCB(TBCA);

	time_t t_filter_start = time(0);
	kappa_filter::run_filter(
		query_kappa_codeseqs2, query_lengths2, nquery, db_ss);

	vector<uint> dbidxs;
	unordered_map<uint, vector<uint> > dbidx_to_qidxs;
	kappa_filter::m_RSB.GetTargetInfo(dbidxs, dbidx_to_qidxs);
	time_t t_filter_end = time(0);
	uint filter_secs = uint(t_filter_end - t_filter_start);

	ProgressLog("Kappa filter %u secs\n", filter_secs);
	}
