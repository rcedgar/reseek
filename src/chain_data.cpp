#include "myutils.h"
#include "chain_data.h"
#include "chaq.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "flat_alphas.h"
#include "flat_bench.h"
#include "flat_nu_aligner.h"
#include "scratch_mem.h"
#include "fan.h"

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values);

static uint32_t s_members;
static const vector<flat_chain_t *> *s_chains;
static chain_data **s_cdvec;
static uint s_nchain;
static atomic<uint> s_next;
static size_t s_memory_bytes;
static size_t s_memory_bytes_per_pos;
static uint8_t *s_memory_base = nullptr;
static size_t s_memory_bytes_total = 0;
static std::atomic<size_t> s_memory_next{0}; // bump offset in bytes

const uint32_t chain_data::m_maxL = 4000;

static void thread_body(uint threadidx)
	{
	const vector<flat_chain_t *> &chains = *s_chains;

	uint scratch_bytes_per_pos =
		4*sizeof(sid_t) +		// chaq_vecs::Xensids
		4*sizeof(uint16_t) +	// chaq_vecs::Xens
		sizeof(uint8_t);		// chaq_vecs::sec32_codeseq
	uint scratch_bytes = chain_data::m_maxL*scratch_bytes_per_pos;
	scratch_mem scratch(scratch_bytes);

	for (;;)
		{
		const uint chainidx = s_next++;
		if (chainidx >= s_nchain)
			return;

		const uint32_t L = chains[chainidx]->get_length(); // or ->m_L if that's your API
		const size_t nbytes = size_t(s_memory_bytes_per_pos) * size_t(L);

		const size_t off = s_memory_next.fetch_add(nbytes, std::memory_order_relaxed);
		asserta(off + nbytes <= s_memory_bytes_total);

		uint8_t *const mem = s_memory_base + off;

		if (threadidx == 0)
			ProgressStep(chainidx, s_nchain, "fill_chain_data_vec");

		s_cdvec[chainidx] =
			chain_data::from_chain(
				*chains[chainidx], s_members,
				mem, nbytes, scratch);
		}
	}

static uint estimate_parasail_profile16_n256_main_bytes(int L)
	{
    uint segLen = (L + 15) / 16;
    return 256*segLen*32;
	}

static void copy_rev_mega_prof(
	uint8_t *rev_prof,
	const uint8_t *prof,
	uint32_t L,
	uint32_t nfeat,
	const uint32_t *alpha_sizes)
	{
	for (uint32_t fi = 0; fi < nfeat; ++fi)
		{
		(void) alpha_sizes[fi];
		for (uint32_t i = 0; i < L; ++i)
			{
			const uint32_t rev_i = L - i - 1;
			rev_prof[size_t(fi) * L + rev_i] = prof[size_t(fi) * L + i];
			}
		}
	}

void chain_data::make_mega_prof(
	const flat_chain_t &chain,
	const sid_t *distmx,
	uint8_t *mega_prof,
	size_t bytes,
	scratch_mem &scratch)
	{
	const uint L = chain.get_length();
	asserta(L > 0);

	const uint nfeat = flat_alphas::m_nfeat;
	asserta(nfeat > 0);
	asserta(bytes >= nfeat*L);

	chaq_vecs cv;
	chaq::fill_chaq_vecs(distmx, L, cv, scratch);

#if DEBUG
	memset(mega_prof, 0xff, nfeat*L);
#endif

	for (uint fi = 0; fi < nfeat; ++fi)
		{
		const FAN fan = flat_alphas::m_fans[fi];
		const uint alpha_size = flat_alphas::m_alpha_sizes[fi];

		uint8_t *get_codeseq_scratch_mem = 0;
		size_t get_codeseq_scratch_bytes = 0;

		uint8_t *codeseq = mega_prof + fi*L;
		uint8_t undef_code = chaq::get_undef_code(fan, alpha_size);
		chaq::fast_get_codeseq(&chain, distmx, &cv, fan, alpha_size,
			codeseq, scratch);

#if DEBUG
		for (uint pos = 0; pos < L; ++pos)
			assert(codeseq[pos] < alpha_size);
#endif
		}
	}

size_t chain_data::get_bytes_per_pos(uint32_t bits)
	{
	const uint32_t M = flat_params::m_distmx_bandwidth;
	const uint32_t nfeat = flat_alphas::m_nfeat;

	size_t bytes = 0;
	if ((bit_distmx & bits) != 0)			{ bytes += M*sizeof(sid_t); }
	if ((bit_mega_pssm & bits) != 0)		{ bytes += flat_alphas::m_sum_alpha_sizes*sizeof(float); }
	if ((bit_mega_pssm_rev & bits) != 0)	{ bytes += flat_alphas::m_sum_alpha_sizes*sizeof(float); }
	if ((bit_mega_prof & bits) != 0)		{ bytes += nfeat; }
	if ((bit_mega_prof_rev & bits) != 0)	{ bytes += nfeat; }
	if ((bit_parasail_prof & bits) != 0)	{ bytes += 0; } // parasail calls malloc
	if ((bit_parasail_prof_rev & bits) != 0){ bytes += 0; } // parasail calls malloc
	if ((bit_nu_codeseq & bits) != 0)		{ bytes += 1; }
	if ((bit_nu_codeseq_rev & bits) != 0)	{ bytes += 1; }
	return bytes;
	}

chain_data *chain_data::from_chain(
	const flat_chain_t &chain,
	uint32_t bits,
	uint8_t *memory,
	size_t memory_bytes,
	scratch_mem &scratch)
	{
	uint8_t *memory_ptr = memory;

	const uint32_t L = chain.get_length();
	asserta(L > 0);

	chain_data *cd = new chain_data;
	cd->m_label = chain.m_label;
	cd->m_chain = &chain;
	cd->m_L = L;

	const uint32_t M = flat_params::m_distmx_bandwidth;
	const uint32_t nfeat = flat_alphas::m_nfeat;

	asserta(bits & bit_distmx);
	cd->m_distmx = (sid_t *) memory_ptr;
	memory_ptr += L*M;
	chaq::fill_distmx(chain.m_xyz->m_data, L, cd->m_distmx);

	const bool want_mega_prof = (bits & bit_mega_prof) != 0;
	const bool want_mega_prof_rev = (bits & bit_mega_prof_rev) != 0;
	const bool want_pssm_fwd = (bits & bit_mega_pssm) != 0;
	const bool want_pssm_rev = (bits & bit_mega_pssm_rev) != 0;
	const bool want_nu = (bits & (bit_nu_codeseq | bit_nu_codeseq_rev)) != 0;
	const bool want_parasail = (bits & (bit_parasail_prof | bit_parasail_prof_rev)) != 0;

	asserta(want_mega_prof);
	size_t prof_bytes = L*nfeat;
	cd->m_mega_prof = memory_ptr;
	memory_ptr += prof_bytes;
	make_mega_prof(chain, cd->m_distmx,
		cd->m_mega_prof, prof_bytes, scratch);

	if (want_mega_prof_rev)
		{
		asserta(cd->m_mega_prof != 0);
		const uint32_t *alpha_sizes = flat_alphas::m_alpha_sizes;
		cd->m_mega_prof_rev = memory_ptr;
		memory_ptr += nfeat*L;
		flat_reverse_profile(
			cd->m_mega_prof, L, nfeat, cd->m_mega_prof_rev);
		}

	if (want_pssm_fwd || want_pssm_rev)
		{
		asserta(want_pssm_fwd);
		asserta(cd->m_mega_prof != 0);
		const uint32_t *alpha_sizes = flat_alphas::m_alpha_sizes;
		const uint nr_pssm_floats = L*flat_alphas::m_sum_alpha_sizes;

		cd->m_mega_pssm = (float *) memory_ptr;
		memory_ptr += nr_pssm_floats*sizeof(float);
		fill_flat_pssm(
			cd->m_mega_prof, L, nfeat, alpha_sizes,
			flat_alphas::m_feature_block_offsets,
			flat_alphas::m_weighted_logoddsvec,
			cd->m_mega_pssm);

		if (want_pssm_rev)
			{
			cd->m_mega_pssm_rev = (float *) memory_ptr;
			memory_ptr += nr_pssm_floats*sizeof(float);
			fill_flat_pssm_reversed(
				cd->m_mega_prof, L, nfeat, alpha_sizes,
				flat_alphas::m_feature_block_offsets,
				flat_alphas::m_weighted_logoddsvec,
				cd->m_mega_pssm_rev);
			}
		}

	if (want_nu || want_parasail)
		{
		asserta(cd->m_mega_prof != 0);

		const uint32_t fi_aa20 = flat_alphas::get_fi(FAN_aa, 20);
		const uint32_t fi_pm2 = flat_alphas::get_fi(FAN_pm, 2);
		const uint32_t fi_sec32 = flat_alphas::get_fi(FAN_sec, 32);

		const uint8_t *prof_aa20 = cd->m_mega_prof + L*size_t(fi_aa20);
		const uint8_t *prof_pm2 = cd->m_mega_prof + L*size_t(fi_pm2);
		const uint8_t *prof_sec32 = cd->m_mega_prof + L*size_t(fi_sec32);

		cd->m_codeseq_nu = memory_ptr;
		memory_ptr += L;

		cd->m_codeseq_nu_rev = memory_ptr;
		memory_ptr += L;

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

			cd->m_codeseq_nu[pos] = code_nu;
			cd->m_codeseq_nu_rev[L - pos - 1] = code_nu;
			}
		}

	if (bits & bit_parasail_prof)
		{
		asserta(cd->m_codeseq_nu != 0);
		cd->m_parasail_prof = parasail_profile_create_avx_256_16(
			(const char *) cd->m_codeseq_nu, L, &flat_nu_aligner::m_matrix);
		}

	if (bits & bit_parasail_prof_rev)
		{
		asserta(cd->m_codeseq_nu_rev != 0);
		cd->m_parasail_prof_rev = parasail_profile_create_avx_256_16(
			(const char *) cd->m_codeseq_nu_rev, L, &flat_nu_aligner::m_matrix);
		}

	size_t bump = size_t(memory_ptr - memory);
	if (bump > memory_bytes)
		Die("bump=%u, memory_bytes %u", uint(bump), uint(memory_bytes));
	return cd;
	}

void chain_data::fill_chain_data_vec(
	const vector<flat_chain_t *> &chains,
	uint32_t bits,
	chain_data **cdvec)
	{
	s_chains = &chains;
	s_nchain = uint(chains.size());
	s_cdvec = cdvec;
	s_members = bits;

	size_t total_length = 0;
	for (auto chain : chains) total_length += chain->m_L;

	s_memory_bytes_per_pos = get_bytes_per_pos(bits);
	s_memory_bytes_total = s_memory_bytes_per_pos*total_length;
	s_memory_base = myalloc64(uint8_t, s_memory_bytes_total);

	const uint nthread = GetRequestedThreadCount();
	ProgressStep(0, s_nchain, "fill_chain_data_vec");

	vector<thread *> ts;
	uint ThreadCount = GetRequestedThreadCount();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(thread_body, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

	ProgressStep(s_nchain-1, s_nchain, "fill_chain_data_vec");
	}

void chain_data::log_mem_stats(chain_data **cdvec, uint nchain)
	{
	uint n_distmx = 0;
	uint n_codeseq_nu = 0;
	uint n_mega_prof = 0;
	uint n_mega_prof_rev = 0;
	uint n_codeseq_nu_rev = 0;
	uint n_mega_pssm = 0;
	uint n_mega_pssm_rev = 0;
	uint n_parasail_prof = 0;
	uint n_parasail_prof_rev = 0;

	size_t bytes_distmx = 0;
	size_t bytes_codeseq_nu = 0;
	size_t bytes_mega_prof = 0;
	size_t bytes_mega_prof_rev = 0;
	size_t bytes_codeseq_nu_rev = 0;
	size_t bytes_mega_pssm = 0;
	size_t bytes_mega_pssm_rev = 0;
	size_t bytes_parasail_prof = 0;
	size_t bytes_parasail_prof_rev = 0;

	const uint32_t M = flat_params::m_distmx_bandwidth;
	const uint32_t nfeat = flat_alphas::m_nfeat;

	for (uint idx = 0; idx < nchain; ++idx)
		{
		const chain_data *cd = cdvec[idx];
		uint L = cd->m_L;
		if (cd->m_distmx != 0)
			{
			++n_distmx;
			bytes_distmx += L*M*sizeof(cd->m_distmx[0]);
			}

		if (cd->m_mega_prof != 0)
			{
			++n_mega_prof;
			bytes_mega_prof += L*nfeat*sizeof(cd->m_mega_prof[0]);
			}

		if (cd->m_mega_prof_rev != 0)
			{
			++n_mega_prof_rev;
			bytes_mega_prof_rev += L*nfeat*sizeof(cd->m_mega_prof[0]);
			}

		if (cd->m_mega_pssm != 0)
			{
			++n_mega_pssm;
			bytes_mega_pssm += L*flat_alphas::m_sum_alpha_sizes*sizeof(m_mega_pssm[0]);
			}

		if (cd->m_mega_pssm_rev != 0)
			{
			++n_mega_pssm_rev;
			bytes_mega_pssm_rev += L*flat_alphas::m_sum_alpha_sizes*sizeof(m_mega_pssm[0]);
			}

		if (cd->m_codeseq_nu != 0)
			{
			++n_codeseq_nu;
			bytes_codeseq_nu += L*sizeof(cd->m_codeseq_nu[0]);
			}

		if (cd->m_codeseq_nu_rev != 0)
			{
			++n_codeseq_nu_rev;
			bytes_codeseq_nu_rev += L*sizeof(cd->m_codeseq_nu_rev[0]);
			}

		if (cd->m_parasail_prof != 0)
			{
			++n_parasail_prof;
			bytes_parasail_prof += estimate_parasail_profile16_n256_main_bytes(L);
			}

		if (cd->m_parasail_prof_rev != 0)
			{
			++n_parasail_prof_rev;
			bytes_parasail_prof_rev += estimate_parasail_profile16_n256_main_bytes(L);
			}
		}
	
	size_t bytes_total = 0;
	bytes_total += bytes_distmx;
	bytes_total += bytes_mega_prof;
	bytes_total += bytes_mega_prof_rev;
	bytes_total += bytes_mega_pssm;
	bytes_total += bytes_mega_pssm_rev;
	bytes_total += bytes_codeseq_nu;
	bytes_total += bytes_codeseq_nu_rev;
	bytes_total += bytes_parasail_prof;
	bytes_total += bytes_parasail_prof_rev;

	Log("%u chains\n", nchain);

#define do(x)	Log("%10u  %12.12s  %6.1f%%  %s\n", \
	n_##x, MemBytesToStr(bytes_##x), GetPct(double(bytes_##x), double(bytes_total)), #x)

	do(distmx);
	do(codeseq_nu);
	do(codeseq_nu_rev);
	do(mega_prof); 
	do(mega_prof_rev);
	do(mega_pssm);
	do(mega_pssm_rev);
	do(parasail_prof);
	do(parasail_prof_rev);
#undef x
	Log("%10u  %12.12s   100.0%%\n",
		nchain, MemBytesToStr(bytes_total));
	}

void cmd_test_chain_data()
	{
	asserta(optset_alphadir);
	asserta(optset_varstr);
	vector<string> param_names;
	vector<float> param_values;

	ParseVarStr(opt(varstr), param_names, param_values);

	vector<string> alpha_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		alpha_names, weights,
		scalar_names, scalar_values);

	const string &alphadir = opt(alphadir);
	flat_alphas::init_from_alphadir(alphadir, alpha_names);
	Paralign::set_final_nu();

	vector<flat_chain_t *> chains;
	read_flat_chains(g_Arg1, chains);
	const uint nchain = uint(chains.size());
	chain_data **cdvec = myalloc(chain_data *, nchain);
	chain_data::fill_chain_data_vec(chains, bits_query, cdvec);
	chain_data::log_mem_stats(cdvec, nchain);
	}
