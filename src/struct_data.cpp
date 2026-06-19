#include "myutils.h"
#include "bcadata.h"
#include "flat_helpers.h"
#include "struct_data.h"
#include "flat_nu_aligner.h"

#include "flat_params.h"

struct_data *BCAData::get_struct_data(
	const flat_params &params,
	uint idx,
	chaq_vecs2 *cv,
	uint8_t *scratch_buffer,
	uint scratch_buffer_bytes) const
	{
	struct_data *sd = new struct_data;

	flat_chain_t *chain = read_flat_chain(idx);
	const uint32_t L = chain->get_length();
	asserta(L > 0);
	asserta(L <= flat_params::m_maxL); // TODO=maxL

	const uint32_t M = params.m_distmx_bandwidth;
	const uint32_t nfeat = params.m_nfeat;

	sid_t *distmx = myalloc(sid_t, L*M);
	chaq::fill_distmx(chain->m_xyz->m_data, L, distmx);

	uint8_t *mega_prof = chaq::make_mega_prof(params, chain, distmx,
		cv, scratch_buffer, scratch_buffer_bytes);

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

	sd->m_chain = chain;
	sd->m_distmx = distmx;
	sd->m_codeseq_nu = codeseq_nu;
	sd->m_codeseq_nu_rev = codeseq_nu_rev;
	sd->m_codeseq_kappa = codeseq_kappa;
	sd->m_mega_prof = mega_prof;
	sd->m_mega_pssm = mega_pssm;
	sd->m_mega_pssm_rev = mega_pssm_rev;
	sd->m_parasail_prof = parasail_prof;
	sd->m_parasail_prof_rev = parasail_prof_rev;
	return sd;
	}

struct_data **BCAData::get_struct_data_vec(const flat_params &params)
	{
	uint nchain = GetChainCount();
	Progress("get_query_data_vec()...");
	uint scratch_buffer_bytes = 2*flat_params::m_maxL;
	uint8_t *scratch_buffer = myalloc(uint8_t, scratch_buffer_bytes);
	chaq_vecs2 cv;
	chaq::alloc_chaq_vecs2(cv, flat_params::m_maxL);
	struct_data **vec = myalloc(struct_data *, nchain);
	for (uint idx = 0; idx < nchain; ++idx)
		vec[idx] = get_struct_data(params, idx,
			&cv, scratch_buffer, scratch_buffer_bytes);
	Progress(" done\n");
	chaq::free_chaq_vecs2(cv);
	myfree(scratch_buffer);
	return vec;
	}

void struct_data::free_struct_data(struct_data *sd)
	{
	delete sd->m_chain;
	myfree(sd->m_distmx);
	myfree(sd->m_codeseq_nu);
	myfree(sd->m_codeseq_nu_rev);
	myfree(sd->m_codeseq_kappa);
	myfree(sd->m_mega_pssm);
	myfree(sd->m_mega_prof);
	myfree(sd->m_mega_pssm_rev);
	parasail_profile_free(sd->m_parasail_prof);
	parasail_profile_free(sd->m_parasail_prof_rev);

	sd->m_chain = 0;
	sd->m_distmx = 0;
	sd->m_codeseq_nu = 0;
	sd->m_codeseq_kappa = 0;
	sd->m_codeseq_nu_rev = 0;
	sd->m_mega_prof = 0;
	sd->m_mega_pssm = 0;
	sd->m_mega_pssm_rev = 0;
	sd->m_parasail_prof = 0;
	sd->m_parasail_prof_rev = 0;
	}
