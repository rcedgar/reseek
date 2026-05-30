#include "myutils.h"
#include "bcadata.h"
#include "flat_helpers.h"
#include "query_data.h"
#include "flat_nu_aligner.h"

query_data *BCAData::get_query_data(
	const flat_params &params,
	uint idx)
	{
	query_data *qd = new query_data;

	flat_chain_t *chain = read_flat_chain(idx);
	const uint32_t L = chain->get_length();
	asserta(L > 0);
	asserta(L <= m_maxL);

	const uint32_t M = params.m_distmx_bandwidth;
	const uint32_t nfeat = params.m_nfeat;

	sid_t *distmx = myalloc(sid_t, L*M);
	chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);

	uint8_t *mega_prof = chaq::make_mega_prof(params, chain, distmx,
		&m_cv, m_scratch_buffer, m_scratch_buffer_bytes);

	const uint32_t *alpha_sizes = params.m_alpha_sizes;
	const uint nr_pssm_floats = L*params.m_sum_alpha_sizes;

	float *mega_pssm = myalloc(float, nr_pssm_floats);
	fill_flat_pssm(
		mega_prof, L, nfeat, alpha_sizes,
		params.m_feature_block_offsets,
		params.m_weighted_logoddsvec,
		mega_pssm);

	float *mega_pssm_rev = myalloc(float, nr_pssm_floats);
	fill_flat_pssm_reversed(
		mega_prof, L, nfeat, alpha_sizes,
		params.m_feature_block_offsets,
		params.m_weighted_logoddsvec,
		mega_pssm_rev);

	const uint32_t fi_aa20 = params.get_fi(FAN_aa, 20);
	const uint32_t fi_pm2 = params.get_fi(FAN_pm, 2);
	const uint32_t fi_sec32 = params.get_fi(FAN_sec, 32);

	const uint8_t *prof_aa20 = mega_prof + L*size_t(fi_aa20);
	const uint8_t *prof_pm2 = mega_prof + L*size_t(fi_pm2);
	const uint8_t *prof_sec32 = mega_prof + L*size_t(fi_sec32);

	uint8_t *codeseq_kappa = myalloc(uint8_t, L);
	uint8_t *codeseq_nu = myalloc(uint8_t, L);
	uint8_t *codeseq_nu_rev = myalloc(uint8_t, L);

	for (uint32_t pos = 0; pos < L; ++pos)
		{
		const uint8_t code_aa20 = prof_aa20[pos];
		const uint8_t code_pm2 = prof_pm2[pos];
		const uint8_t code_sec32 = prof_sec32[pos];

		assert(code_aa20 < 20);
		assert(code_pm2 < 2);
		assert(code_sec32 < 32);

		const uint8_t code_aa4 = chaq::m_aacode2aa4code[code_aa20];
		const uint8_t code_nu = uint8_t(code_aa4 + 4*code_pm2 + 4*2*code_sec32);
		assert(code_nu < 256);

		codeseq_nu[pos] = code_nu;
		codeseq_nu_rev[L-pos-1] = code_nu;
		
		chaq::codeseq_nu_to_kappa(codeseq_nu, L, codeseq_kappa, L);
		}

	parasail_profile_t *parasail_prof = parasail_profile_create_avx_256_16(
		(const char *) codeseq_nu, L, &flat_nu_aligner::m_matrix);

	parasail_profile_t *parasail_prof_rev = parasail_profile_create_avx_256_16(
		(const char *) codeseq_nu_rev, L, &flat_nu_aligner::m_matrix);

	qd->m_chain = chain;
	qd->m_distmx = distmx;
	qd->m_codeseq_nu = codeseq_nu;
	qd->m_codeseq_kappa = codeseq_kappa;
	qd->m_mega_prof = mega_prof;
	qd->m_mega_pssm = mega_pssm;
	qd->m_mega_pssm_rev = mega_pssm_rev;
	qd->m_parasail_prof = parasail_prof;
	qd->m_parasail_prof_rev = parasail_prof_rev;
	return qd;
	}

query_data **BCAData::get_query_data_vec(const flat_params &params)
	{
	uint nchain = GetChainCount();
	Progress("get_query_data_vec()...");
	query_data **vec = myalloc(query_data *, nchain);
	for (uint idx = 0; idx < nchain; ++idx)
		vec[idx] = get_query_data(params, idx);
	Progress(" done\n");
	return vec;
	}
