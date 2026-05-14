#include "myutils.h"
#include "flat_bench2.h"
#include "thread_affinity.h"
#include "flat_helpers.h"
#include "flat_alphas.h"
#include "paralign.h"

chain_data **flat_bench2::m_cdvec;
uint flat_bench2::m_maxL = 4000;

void flat_bench2::search(uint nthread, bool pin_threads)
	{
	FastBench::Alloc();

	const uint NQ = SIZE(m_Labels);
	const uint PairCount = triangle_get_K(NQ);

	ProgressStep(0, PairCount, "Aligning");
	m_next_pairidx = 0;
	m_aligned_pair_count = 0;

	thread_affinity ta;
	vector<thread *> ts;
	for (uint threadidx = 0; threadidx < nthread; ++threadidx)
		{
		thread *t = new thread(static_thread_body, this, threadidx);
		if (pin_threads)
			ta.pinThread(*t, threadidx);
		ts.push_back(t);
		}
	for (uint threadidx = 0; threadidx < nthread; ++threadidx)
		ts[threadidx]->join();
	for (uint threadidx = 0; threadidx < nthread; ++threadidx)
		delete ts[threadidx];
	ProgressStep(PairCount-1, PairCount, "Aligning");
	}

void flat_bench2::static_thread_body(flat_bench2 *FB, uint threadidx)
	{
	FB->thread_body(threadidx);
	}

void flat_bench2::thread_body(uint threadidx)
	{
	const uint NQ = SIZE(m_Labels);
	const uint npair = triangle_get_K(NQ);

	for (;;)
		{
		uint pairidx = m_next_pairidx++;
		if (pairidx >= npair)
			return;
		uint progress_count = m_aligned_pair_count++;
		if (threadidx == 0 && progress_count%1000 == 0)
			ProgressStep(progress_count, npair, "Aligning");
		align_pair(pairidx);
		}
	}

void flat_bench2::load_chains(const vector<flat_chain_t *> &chains)
	{
	uint nchain = uint(chains.size());
	m_cdvec = myalloc(chain_data *, nchain);
	chain_data::fill_chain_data_vec(chains, bits_query, m_cdvec);
	}

void flat_bench2::align_pair(uint pairidx)
	{
	uint NQ = uint(m_Labels.size());
	uint i, j;
	triangle_k_to_ij(pairidx, NQ, i, j);

	const chain_data *cd_i = m_cdvec[i];
	const chain_data *cd_j = m_cdvec[j];

	const uint L_i = cd_i->m_L;
	const uint L_j = cd_j->m_L;
	asserta(L_i <= m_maxL);
	asserta(L_j <= m_maxL);

	const uint8_t *prof_i = cd_i->m_mega_prof;
	const float *pssm_j = cd_j->m_mega_pssm;

	uint lo_i, lo_j, ncol;
	float score = sw_flat_pssm(
		m_scratch_rows, m_TB, m_scratch_pssms,
		prof_i, L_i,
		pssm_j, L_j, 
		flat_alphas::m_feature_block_offsets,
		flat_alphas::m_nfeat,
		-flat_params::m_open, 
		-flat_params::m_ext,
		lo_i, lo_j, m_path_buffer, ncol);

	m_Scores[pairidx] = score;
	}

void flat_bench2::update_params(
	const vector<string> &names,
	const vector<float> &values)
	{
	vector<string> alpha_names;
	vector<float> weights;
	vector<string> scalar_names;
	vector<float> scalar_values;
	flat_classify_params(
		names, values, alpha_names,
		weights, scalar_names, scalar_values);

	flat_params::set_params(scalar_names, scalar_values);

	uint n = SIZE(alpha_names);
	asserta(SIZE(weights) == n);
	unordered_map<string, float> NameToWeight;
	for (uint i = 0; i < n; ++i)
		{
		const string &name = alpha_names[i];
		if (NameToWeight.find(name) != NameToWeight.end())
			Die("Dupe name in spec '%s'", name.c_str());
		NameToWeight[name] = weights[i];
		}
	flat_alphas::apply_weights(NameToWeight);
	if (flat_params::need_self())
		Warning("self scores not implemented");
	}

void cmd_flat_bench2()
	{
	Paralign::set_final_nu();

	vector<string> param_names;
	vector<float> param_values;
	parse_varstr(opt(varstr), param_names, param_values);

	vector<string> alpha_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_classify_params(
		param_names, param_values,
		alpha_names, weights,
		scalar_names, scalar_values);

	const string &alphadir = opt(alphadir);
	flat_alphas::init_from_alphadir(alphadir, alpha_names);

	vector<flat_chain_t *> chains;
	read_flat_chains(g_Arg1, chains);

	flat_bench2 FB;
	FB.ReadLookup(opt(lookup));
	FB.update_params(param_names, param_values);
	FB.load_chains(chains);

	uint nthread = GetRequestedThreadCount();
	thread_affinity ta;
	bool pin = opt(no_thread_pin) ? false : ta.shouldPin(nthread);
	FB.search(nthread, pin);
	FB.SetScoreOrder();
	FB.Bench();
	FB.WriteHits(opt(output), opt(include_self), opt(triangle));
	}
