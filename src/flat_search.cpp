#include "myutils.h"
#include "kappa_filter.h"
#include "kappa_mermx.h"
#include "kappa_dex.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "nu_filter.h"

#define WRITE_TARGET_INFO	0

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
	struct_data **struct_data_vec = QBCA.get_struct_data_vec(params);

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
	uint8_t **query_mega_profs = myalloc(uint8_t *, nquery);
	const float **query_mega_pssms = myalloc(const float *, nquery);
	const float **query_mega_pssm_revs = myalloc(const float *, nquery);
	uint *query_lengths = myalloc(uint, nquery);
	const sid_t **query_distmxs = myalloc(const sid_t *, nquery);
	parasail_profile_t **query_parasail_profs = myalloc(parasail_profile_t *, nquery);
	parasail_profile_t **query_parasail_prof_revs = myalloc(parasail_profile_t *, nquery);

	for (uint chainidx = 0; chainidx < nquery; ++chainidx)
		{
		uint L = struct_data_vec[chainidx]->m_chain->get_length();
		query_lengths[chainidx] = L;
		query_nu_codeseqs[chainidx] = struct_data_vec[chainidx]->m_codeseq_nu;
		query_kappa_codeseqs[chainidx] = struct_data_vec[chainidx]->m_codeseq_kappa;
		query_parasail_profs[chainidx] = struct_data_vec[chainidx]->m_parasail_prof;
		query_parasail_prof_revs[chainidx] = struct_data_vec[chainidx]->m_parasail_prof_rev;
		query_mega_profs[chainidx] = struct_data_vec[chainidx]->m_mega_prof;
		query_mega_pssms[chainidx] = struct_data_vec[chainidx]->m_mega_pssm;
		query_mega_pssm_revs[chainidx] = struct_data_vec[chainidx]->m_mega_pssm_rev;
		query_distmxs[chainidx] = struct_data_vec[chainidx]->m_distmx;
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
#if WRITE_TARGET_INFO
	{
	FILE *f = CreateStdioFile("flat_search_target_info.tmp");
	for (size_t i = 0; i < dbidxs.size(); ++i)
		{
		uint dbidx = dbidxs[i];
		const vector<uint> &qidxs = dbidx_to_qidxs[dbidx];
		for (auto qidx : qidxs)
			{
			const char *qlabel = QBCA.m_Labels[qidx].c_str();
			const char *dblabel = DBBCA.m_Labels[dbidx].c_str();
			fprintf(f, "%s\t%s\n", qlabel, dblabel);
			}
		}
	}
#endif
	time_t t_kappa_filter_end = time(0);
	uint kappa_filter_secs = uint(t_kappa_filter_end - t_kappa_filter_start);
	ProgressLog("Kappa filter %u secs\n", kappa_filter_secs);

	nu_filter::set_params(params);
	nu_filter::set_query_data(
		QBCA.m_Labels,
		query_parasail_profs,
		query_parasail_prof_revs,
		query_mega_pssms,
		query_mega_pssm_revs,
		query_distmxs,
		query_lengths,
		nquery);

	nu_filter::set_query_self_rev_scores(query_nu_codeseqs);
	nu_filter::set_query_mega_self_rev_scores(query_mega_profs);
	nu_filter::run_filter(DBBCA, dbidxs, dbidx_to_qidxs);
	time_t t_nu_filter_end = time(0);
	uint nu_filter_secs = uint(t_nu_filter_end - t_kappa_filter_end);
	ProgressLog("Nu filter %u secs\n", kappa_filter_secs);
	}
