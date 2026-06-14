#include "myutils.h"
#include "kappa_filter.h"
#include "kappa_mermx.h"
#include "kappa_dex.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "reseeker.h"

void cmd_flat_search_all()
	{
	const string &QFN = g_Arg1;
	const string &DBFN = opt(db);

	BCAData QBCA;
	BCAData DBBCA;

	QBCA.Open(QFN);
	asserta(QBCA.m_HasNuSequences);

	flat_params params_fold;
	flat_params params_sf;
	flat_params params_fam;
	flat_init_reseeker_params(params_fold, params_sf, params_fam);
	params_fold.logme();
	struct_data **struct_data_vec = QBCA.get_struct_data_vec(params_fold);

	DBBCA.Open(DBFN);

	const uint nquery = QBCA.GetChainCount();
	const uint TSeqCount = QBCA.GetChainCount();

	const flat_chain_t **query_chains = myalloc(const flat_chain_t *, nquery);
	uint8_t **query_nu_codeseqs = myalloc(uint8_t *, nquery);
	uint8_t **query_mega_profs = myalloc(uint8_t *, nquery);
	const float **query_mega_pssms_fold = myalloc(const float *, nquery);
	const float **query_mega_pssm_revs_fold = myalloc(const float *, nquery);
	const float **query_mega_pssms_sf = 0;
	const float **query_mega_pssm_revs_sf = 0;
	const float **query_mega_pssms_fam = 0;
	const float **query_mega_pssm_revs_fam = 0;
	uint *query_lengths = myalloc(uint, nquery);
	const sid_t **query_distmxs = myalloc(const sid_t *, nquery);
	parasail_profile_t **query_parasail_profs = myalloc(parasail_profile_t *, nquery);
	parasail_profile_t **query_parasail_prof_revs = myalloc(parasail_profile_t *, nquery);

	for (uint chainidx = 0; chainidx < nquery; ++chainidx)
		{
		uint L = struct_data_vec[chainidx]->m_chain->get_length();
		query_chains[chainidx] = struct_data_vec[chainidx]->m_chain;
		query_lengths[chainidx] = L;
		query_nu_codeseqs[chainidx] = struct_data_vec[chainidx]->m_codeseq_nu;
		query_parasail_profs[chainidx] = struct_data_vec[chainidx]->m_parasail_prof;
		query_parasail_prof_revs[chainidx] = struct_data_vec[chainidx]->m_parasail_prof_rev;
		query_mega_profs[chainidx] = struct_data_vec[chainidx]->m_mega_prof;
		query_mega_pssms_fold[chainidx] = struct_data_vec[chainidx]->m_mega_pssm;
		query_mega_pssm_revs_fold[chainidx] = struct_data_vec[chainidx]->m_mega_pssm_rev;
		query_distmxs[chainidx] = struct_data_vec[chainidx]->m_distmx;
		}

	flat_build_query_mega_pssms(
		nquery, query_mega_profs, query_lengths,
		params_sf, query_mega_pssms_sf, query_mega_pssm_revs_sf);
	flat_build_query_mega_pssms(
		nquery, query_mega_profs, query_lengths,
		params_fam, query_mega_pssms_fam, query_mega_pssm_revs_fam);

	reseeker::set_params_fold(&params_fold);
	reseeker::set_params_sf(&params_sf);
	reseeker::set_params_fam(&params_fam);
	reseeker::set_query_data(
		query_chains,
		QBCA.m_Labels,
		query_parasail_profs,
		query_parasail_prof_revs,
		query_mega_pssms_fold,
		query_mega_pssm_revs_fold,
		query_mega_pssms_sf,
		query_mega_pssm_revs_sf,
		query_mega_pssms_fam,
		query_mega_pssm_revs_fam,
		query_distmxs,
		query_lengths,
		nquery);

	reseeker::set_query_self_rev_scores(query_nu_codeseqs);
	reseeker::set_query_mega_self_rev_scores(query_mega_profs);
	reseeker::search_all_vs_all(DBBCA);
	}
