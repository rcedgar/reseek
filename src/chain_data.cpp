#include "myutils.h"
#include "chain_data.h"
#include "chaq.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "flat_nu_aligner.h"
#include "scratch_mem.h"
#include "fan.h"

static uint32_t s_bits;
static const vector<flat_chain_t *> *s_chains;
static chain_data **s_cdvec;
static uint s_nchain;
static atomic<uint> s_next;
static size_t s_mem_bytes;
static uint8_t *s_mem_base = nullptr;
static size_t s_mem_bytes_total = 0;
static atomic<size_t> s_mem_next{0}; // bump offset in bytes
static const flat_params *s_params = 0;

static size_t s_mem_bytes_per_pos;
static size_t s_scratch_bytes_per_pos;

static void thread_body(uint threadidx)
	{
	const vector<flat_chain_t *> &chains = *s_chains;

	scratch_mem scratch(flat_params::m_maxL*s_scratch_bytes_per_pos);

	uint scratch_buffer_bytes = 2*flat_params::m_maxL;
	uint8_t *scratch_buffer = myalloc(uint8_t, scratch_buffer_bytes);
	chaq_vecs2 cv;
	chaq::alloc_chaq_vecs2(cv, flat_params::m_maxL);
	for (;;)
		{
		const uint chainidx = s_next++;
		if (chainidx >= s_nchain)
			return;

		const uint32_t L = chains[chainidx]->get_length();
		asserta(L <= flat_params::m_maxL); // TODO=maxL
		const size_t nbytes = L*s_mem_bytes_per_pos;

		const size_t off = s_mem_next.fetch_add(nbytes, std::memory_order_relaxed);
		asserta(off + nbytes <= s_mem_bytes_total);

		uint8_t *ptr_mem = s_mem_base + off;
		scratch_mem mem(ptr_mem, nbytes);

		if (threadidx == 0 && chainidx + 1 < s_nchain)
			ProgressStep(chainidx, s_nchain, "fill_chain_data_vec");

		s_cdvec[chainidx] =
			chain_data::from_chain(
				*s_params,
				*chains[chainidx], s_bits,
				mem, scratch_buffer, scratch_buffer_bytes, &cv);

		scratch.reset();
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
	const flat_params &params,
	const flat_chain_t &chain,
	const sid_t *distmx,
	uint8_t *mega_prof,
	uint mega_prof_bytes,
	uint8_t *scratch_buffer,
	uint scratch_buffer_bytes,
	chaq_vecs2 *cv)
	{
	const uint L = chain.get_length();
	asserta(L > 0);

	const uint nfeat = params.m_nfeat;
	asserta(nfeat > 0);
	asserta(mega_prof_bytes >= nfeat*L);

	chaq::fill_chaq_vecs2(distmx, L, *cv);

#if DEBUG
	memset(mega_prof, 0xff, nfeat*L);
#endif

	for (uint fi = 0; fi < nfeat; ++fi)
		{
		const FAN fan = params.m_fans[fi];
		const uint alpha_size = params.m_alpha_sizes[fi];

		uint8_t *codeseq = mega_prof + fi*L;
		uint8_t undef_code = chaq::get_undef_code(fan, alpha_size);
		chaq::fast_get_codeseq(
			params, &chain, distmx,
			cv, fan, alpha_size, codeseq,
			scratch_buffer, scratch_buffer_bytes);

#if DEBUG
		for (uint pos = 0; pos < L; ++pos)
			assert(codeseq[pos] < alpha_size);
#endif
		}
	}

void chain_data::get_from_chain_bytes_per_pos(
	const flat_params &params,
	uint32_t bits,
	size_t &mem_bytes_per_pos,
	size_t &scratch_bytes_per_pos)
	{
	const uint32_t M = params.m_distmx_bandwidth;
	const uint32_t nfeat = params.m_nfeat;

	size_t fill_chaq_vecs_bytes_per_pos =
		chaq::get_fill_chaq_vecs_bytes_per_pos();

	size_t fast_get_codeseq_scratch_bytes_per_pos =
		chaq::get_fast_get_codeseq_scratch_bytes_per_pos();

	uint n_distmx, n_mega_prof, n_mega_pssm, n_nu_codeseq;
	get_object_counts(bits, n_distmx, n_mega_prof, n_mega_pssm, n_nu_codeseq);

	const uint nr_pssm_floats_per_pos = params.m_sum_alpha_sizes;

	mem_bytes_per_pos = 0;
	mem_bytes_per_pos += n_distmx*M*sizeof(sid_t);
	mem_bytes_per_pos += n_mega_prof*nfeat;
	mem_bytes_per_pos += n_mega_pssm*nr_pssm_floats_per_pos*sizeof(float);
	mem_bytes_per_pos += n_nu_codeseq;

	scratch_bytes_per_pos = 0;
	scratch_bytes_per_pos += fill_chaq_vecs_bytes_per_pos;

	if (n_nu_codeseq > 0) scratch_bytes_per_pos += fast_get_codeseq_scratch_bytes_per_pos;
	}

void chain_data::update_pssms(const flat_params &params, chain_data **cdvec, uint n)
	{
	asserta(cdvec != 0);
	for (uint i = 0; i < n; ++i) update_pssms_cd(params, cdvec[i]);
	}

void chain_data::update_pssms_cd(const flat_params &params, chain_data *cd)
	{
	asserta(cd != 0);
	const uint L = cd->m_L;
	const uint32_t *alpha_sizes = params.m_alpha_sizes;
	const uint nr_pssm_floats = L*params.m_sum_alpha_sizes;
	const uint nfeat = params.m_nfeat;

	fill_flat_pssm(
		cd->m_mega_prof, L, nfeat, alpha_sizes,
		params.m_feature_block_offsets,
		params.m_weighted_logoddsvec,
		cd->m_mega_pssm);

	if (cd->m_mega_pssm_rev)
		fill_flat_pssm_reversed(
			cd->m_mega_prof, L, nfeat, alpha_sizes,
			params.m_feature_block_offsets,
			params.m_weighted_logoddsvec,
			cd->m_mega_pssm_rev);
	}

chain_data *chain_data::from_chain(
	const flat_params &params,
	const flat_chain_t &chain,
	uint32_t bits,
	scratch_mem &mem,
	uint8_t *scratch_buffer,
	uint scratch_buffer_bytes,
	chaq_vecs2 *cv)
	{
	asserta(bits & bit_distmx);

	const uint32_t L = chain.get_length();
	asserta(L > 0);
	asserta(L <= flat_params::m_maxL); // TODO=maxL

	chain_data *cd = new chain_data;
	cd->m_label = chain.m_label;
	cd->m_chain = &chain;
	cd->m_L = L;

	const uint32_t M = params.m_distmx_bandwidth;
	const uint32_t nfeat = params.m_nfeat;

	asserta(bits & bit_distmx);
	cd->m_distmx = mem.get<sid_t>(L*M);
	chaq::fill_distmx(chain.m_xyz->m_data, L, cd->m_distmx);

	const bool want_mega_prof = (bits & bit_mega_prof) != 0;
	const bool want_pssm_fwd = (bits & bit_mega_pssm) != 0;
	const bool want_pssm_rev = (bits & bit_mega_pssm_rev) != 0;
	const bool want_nu_codeseq = (bits & bit_nu_codeseq) != 0;

	asserta(want_mega_prof);
	uint prof_bytes = L*nfeat;
	cd->m_mega_prof = mem.get<uint8_t>(L*nfeat);
	make_mega_prof(params, chain, cd->m_distmx,
		cd->m_mega_prof, prof_bytes,
		scratch_buffer, scratch_buffer_bytes, cv);

	if (want_pssm_fwd || want_pssm_rev)
		{
		asserta(want_pssm_fwd);
		asserta(cd->m_mega_prof != 0);
		const uint32_t *alpha_sizes = params.m_alpha_sizes;
		const uint nr_pssm_floats = L*params.m_sum_alpha_sizes;

		cd->m_mega_pssm = mem.get<float>(nr_pssm_floats);
		fill_flat_pssm(
			cd->m_mega_prof, L, nfeat, alpha_sizes,
			params.m_feature_block_offsets,
			params.m_weighted_logoddsvec,
			cd->m_mega_pssm);

		if (want_pssm_rev)
			{
			cd->m_mega_pssm_rev = mem.get<float>(nr_pssm_floats);
			fill_flat_pssm_reversed(
				cd->m_mega_prof, L, nfeat, alpha_sizes,
				params.m_feature_block_offsets,
				params.m_weighted_logoddsvec,
				cd->m_mega_pssm_rev);
			}
		}

	if (want_nu_codeseq)
		{
		asserta(cd->m_mega_prof != 0);

		const uint32_t fi_aa20 = params.get_fi(FAN_aa, 20);
		const uint32_t fi_pm2 = params.get_fi(FAN_pm, 2);
		const uint32_t fi_sec32 = params.get_fi(FAN_sec, 32);

		const uint8_t *prof_aa20 = cd->m_mega_prof + L*size_t(fi_aa20);
		const uint8_t *prof_pm2 = cd->m_mega_prof + L*size_t(fi_pm2);
		const uint8_t *prof_sec32 = cd->m_mega_prof + L*size_t(fi_sec32);

		cd->m_codeseq_nu = mem.get<uint8_t>(L);
		cd->m_codeseq_nu_rev = mem.get<uint8_t>(L);

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
			cd->m_codeseq_nu_rev[L-pos-1] = code_nu;
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
		asserta(cd->m_codeseq_nu != 0);
		cd->m_parasail_prof_rev = parasail_profile_create_avx_256_16(
			(const char *) cd->m_codeseq_nu_rev, L, &flat_nu_aligner::m_matrix);
		}

	return cd;
	}

void chain_data::get_object_counts(
	uint32_t bits,
	uint &n_distmx,
	uint &n_mega_prof,
	uint &n_mega_pssm,
	uint &n_nu_codeseq)
	{
	n_distmx = 0;
	n_mega_prof = 0;
	n_mega_pssm = 0;
	n_nu_codeseq = 0;

	if (bits & bit_distmx) ++n_distmx;
	if (bits & bit_mega_prof) ++n_mega_prof;
	if (bits & bit_mega_pssm) ++n_mega_pssm;
	if (bits & bit_mega_pssm_rev) ++n_mega_pssm;
	if (bits & bit_nu_codeseq) ++n_nu_codeseq;
	if (bits & bit_nu_codeseq_rev) ++n_nu_codeseq;
	}

void chain_data::fill_chain_data_vec(
	const flat_params &params,
	const vector<flat_chain_t *> &chains,
	uint32_t bits,
	chain_data **cdvec)
	{
	s_chains = &chains;
	s_nchain = uint(chains.size());
	s_cdvec = cdvec;
	s_bits = bits;
	s_params = &params;

	size_t total_length = 0;
	for (auto chain : chains) total_length += chain->m_L;

	get_from_chain_bytes_per_pos(
		params, bits, s_mem_bytes_per_pos, s_scratch_bytes_per_pos);

	s_mem_bytes_total = s_mem_bytes_per_pos*total_length;
	s_mem_base = myalloc64(uint8_t, s_mem_bytes_total);
#if DEBUG
	memset(s_mem_base, 0xff, s_mem_bytes_total);
#endif

	const uint nthread = GetRequestedThreadCount();
	ProgressStep(0, s_nchain, "fill_chain_data_vec");

	s_next = 0;
	s_mem_next = 0;
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

void chain_data::log_mem_stats(
	const flat_params &params,
	chain_data **cdvec,
	uint nchain)
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

	const uint32_t M = params.m_distmx_bandwidth;
	const uint32_t nfeat = params.m_nfeat;

	for (uint idx = 0; idx < nchain; ++idx)
		{
		const chain_data *cd = cdvec[idx];
		uint L = cd->m_L;
		asserta(L <= flat_params::m_maxL); // TODO=maxL
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

		if (cd->m_mega_pssm != 0)
			{
			++n_mega_pssm;
			bytes_mega_pssm += L*params.m_sum_alpha_sizes*sizeof(m_mega_pssm[0]);
			}

		if (cd->m_mega_pssm_rev != 0)
			{
			++n_mega_pssm_rev;
			bytes_mega_pssm_rev += L*params.m_sum_alpha_sizes*sizeof(m_mega_pssm[0]);
			}

		if (cd->m_codeseq_nu != 0)
			{
			++n_codeseq_nu;
			bytes_codeseq_nu += L*sizeof(cd->m_codeseq_nu[0]);
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
	n_##x, MemBytesToStr(double(bytes_##x)), GetPct(double(bytes_##x), double(bytes_total)), #x)

	do(codeseq_nu);
	do(codeseq_nu_rev);
	do(mega_prof); 
	do(mega_prof_rev);
	do(distmx);
	do(parasail_prof);
	do(parasail_prof_rev);
	do(mega_pssm);
	do(mega_pssm_rev);
#undef x
	Log("%10u  %12.12s   100.0%%\n",
		nchain, MemBytesToStr(double(bytes_total)));
	}

void chain_data::write_fastas(
	const flat_params &params,
	const string &fnprefix,
	const chain_data *const *cdvec,
	uint nchain)
	{
	if (fnprefix == "") return;
	const uint nfeat = params.m_nfeat;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		FAN fan = params.m_fans[fi];
		uint alpha_size = params.m_alpha_sizes[fi];
		asserta(alpha_size <= 32);
		string fn;
		Ps(fn, "%s.%s%u.fasta",
			fnprefix.c_str(),
			string(FAN2str(fan)).c_str(),
			alpha_size);
		ProgressStep(fi, nfeat, "%s", fn.c_str());

		FILE *f = CreateStdioFile(fn);
		for (uint chainidx = 0; chainidx < nchain; ++chainidx)
			{
			const chain_data *cd = cdvec[chainidx];
			const string &label = cd->m_label;
			const uint L = cd->m_L;
			const uint8_t *mega_prof = cd->m_mega_prof;
			asserta(mega_prof != 0);
			string seq;
			seq.reserve(L);
			const uint8_t *code2char =
				(fan == FAN_aa && alpha_size == 20 ?
				g_LetterToCharAmino : g_LetterToCharMu);
			for (uint pos = 0; pos < L; ++pos)
				{
				uint8_t code = mega_prof[L*fi + pos];
				asserta(code < alpha_size);
				char c = code2char[code];
				seq += c;
				}
			SeqToFasta(f, label, seq);
			}
		CloseStdioFile(f);
		}
	}

#if 0
void cmd_test_chain_data()
	{
	//asserta(optset_varstr);
	vector<string> param_names;
	vector<float> param_values;

	parse_varstr(opt(varstr), param_names, param_values);

	vector<string> alpha_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		alpha_names, weights,
		scalar_names, scalar_values);

	const string &alphadir = opt(alphadir);
	flat_params params;
	params.init_from_alphadir(alphadir, alpha_names);
	s_params = &params;
	Paralign::set_final_nu();

	vector<flat_chain_t *> chains;
	read_flat_chains(g_Arg1, chains);
	const uint nchain = uint(chains.size());

	chain_data **cdvec = myalloc(chain_data *, nchain);
	chain_data::fill_chain_data_vec(params, chains, bits_query, cdvec);
	chain_data::log_mem_stats(params, cdvec, nchain);
	chain_data::write_fastas(params, opt(fasta), cdvec, nchain);
	}
#endif
