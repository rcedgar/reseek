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
	m_aln_count = 0;
	m_mu_fwd_reject_count = 0;
	m_mu_combined_reject_count = 0;

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

	uint nfeat = flat_alphas::m_nfeat;
	asserta(nfeat > 0);

	flat_bench2_thread_data TD(m_maxL, nfeat);
	for (;;)
		{
		uint pairidx = m_next_pairidx++;
		if (pairidx >= npair)
			return;
		uint progress_count = m_aln_count++;
		if (threadidx == 0 && progress_count%1000 == 0)
			ProgressStep(progress_count, npair, "Aligning");
		align_pair(pairidx, TD);
		}
	}

void flat_bench2::set_nu_self_rev_scores()
	{
	if (!flat_params::need_nu_self())
		{
		asserta(m_nu_self_rev_scores == 0);
		return;
		}

	flat_bench2_thread_data TD(m_maxL, flat_alphas::m_nfeat);
	const uint ndom = m_look->get_ndom();
	if (m_nu_self_rev_scores == 0)
		m_nu_self_rev_scores = myalloc(float, ndom);

	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		ProgressStep(domidx, ndom, "Nu self-rev");
		const chain_data *cd = m_cdvec[domidx];
		const int open = Paralign::m_Open;
		const int ext = Paralign::m_Ext;

		if (TD.m_parasail_result != 0)
			parasail_result_free(TD.m_parasail_result);
		parasail_profile_t *prof = cd->m_parasail_prof;
		asserta(prof != 0);
		const uint8_t *codeseq_nu_rev = cd->m_codeseq_nu_rev;
		asserta(codeseq_nu_rev != 0);
		TD.m_parasail_result = parasail_sw_striped_profile_avx2_256_16(
			prof, (const char *) codeseq_nu_rev, cd->m_L, open, ext);
		asserta(!(TD.m_parasail_result->flag & PARASAIL_FLAG_SATURATED));
		m_nu_self_rev_scores[domidx] = float(TD.m_parasail_result->score);
		}
	}

void flat_bench2::set_self_rev_scores()
	{
	if (!flat_params::need_self())
		{
		asserta(m_self_rev_scores == 0);
		return;
		}

	flat_bench2_thread_data TD(m_maxL, flat_alphas::m_nfeat);
	const uint ndom = m_look->get_ndom();
	if (m_self_rev_scores == 0)
		m_self_rev_scores = myalloc(float, ndom);

	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		ProgressStep(domidx, ndom, "Mega self-rev");
		const chain_data *cd = m_cdvec[domidx];
		const uint8_t *prof = cd->m_mega_prof_rev;
		const float *pssm = cd->m_mega_pssm;
		asserta(prof != 0);
		asserta(pssm != 0);
		const uint L = cd->m_L;
		uint lo_i, lo_j, ncol;
		float score = sw_flat_pssm(
			TD.m_scratch_rows, TD.m_TB, TD.m_scratch_pssms,
			prof, L, pssm, L, 
			flat_alphas::m_feature_block_offsets,
			flat_alphas::m_nfeat,
			-flat_params::m_open, 
			-flat_params::m_ext,
			lo_i, lo_j, TD.m_path_buffer, ncol);
		m_self_rev_scores[domidx] = score;
		}
	}

void flat_bench2::load_chains(const vector<flat_chain_t *> &chains)
	{
	uint nchain = uint(chains.size());
	m_cdvec = myalloc(chain_data *, nchain);
	vector<flat_chain_t *> sorted_chains;
	m_look->sort_chains(chains, sorted_chains);
	chain_data::fill_chain_data_vec(sorted_chains, bits_query, m_cdvec);
	}

void flat_bench2::align_pair(
	uint pairidx, flat_bench2_thread_data &TD)
	{
	uint NQ = uint(m_Labels.size());
	uint i, j;
	triangle_k_to_ij(pairidx, NQ, i, j);

	const chain_data *cd_i = m_cdvec[i];
	const chain_data *cd_j = m_cdvec[j];

	string label_i = cd_i->m_label;
	string label_j = cd_j->m_label;
	trunc_label(label_i);
	trunc_label(label_j);
	asserta(label_i == m_look->get_dom(i));
	asserta(label_j == m_look->get_dom(j));

	const uint L_i = cd_i->m_L;
	const uint L_j = cd_j->m_L;
	asserta(L_i <= m_maxL);
	asserta(L_j <= m_maxL);

	if (flat_params::m_min_nu_fwd_score > 0)
		{
		const int open = Paralign::m_Open;
		const int ext = Paralign::m_Ext;

		if (TD.m_parasail_result != 0)
			parasail_result_free(TD.m_parasail_result);
		parasail_profile_t *prof_i = cd_i->m_parasail_prof;
		asserta(prof_i != 0);
		const uint8_t *codeseq_nu_j = cd_j->m_codeseq_nu;
		TD.m_parasail_result = parasail_sw_striped_profile_avx2_256_16(
			prof_i, (const char *) codeseq_nu_j, L_j, open, ext);
		asserta(!(TD.m_parasail_result->flag & PARASAIL_FLAG_SATURATED));
		int fwd_score = TD.m_parasail_result->score;
		if (fwd_score < flat_params::m_min_nu_fwd_score)
			{
			++m_mu_fwd_reject_count;
			return;
			}

		parasail_result_free(TD.m_parasail_result);
		parasail_profile_t *prof_i_rev = cd_i->m_parasail_prof_rev;
		asserta(prof_i_rev != 0);
		TD.m_parasail_result = parasail_sw_striped_profile_avx2_256_16(
			prof_i, (const char *) codeseq_nu_j, L_j, open, ext);
		asserta(!(TD.m_parasail_result->flag & PARASAIL_FLAG_SATURATED));
		int rev_score = TD.m_parasail_result->score;

		float self_score = (m_nu_self_rev_scores[i] + 
			m_nu_self_rev_scores[j])/2.0f;

		const float revw = flat_params::m_nu_filter_rev_w;
		const float selfw = flat_params::m_nu_filter_self_w;

		float nu_combined_score =
			float(fwd_score) -
			selfw*self_score -
			revw*float(rev_score);
		if (nu_combined_score < flat_params::m_min_nu_combined_score)
			{
			++m_mu_combined_reject_count;
			return;
			}
		}
	
	const uint8_t *prof_i = cd_i->m_mega_prof;
	const float *pssm_j = cd_j->m_mega_pssm;

	uint lo_i, lo_j, ncol;
	float score = sw_flat_pssm(
		TD.m_scratch_rows, TD.m_TB, TD.m_scratch_pssms,
		prof_i, L_i,
		pssm_j, L_j, 
		flat_alphas::m_feature_block_offsets,
		flat_alphas::m_nfeat,
		-flat_params::m_open, 
		-flat_params::m_ext,
		lo_i, lo_j, TD.m_path_buffer, ncol);
	const string path = string(TD.m_path_buffer);

	const sid_t *distmx_i = cd_i->m_distmx;
	const sid_t *distmx_j = cd_j->m_distmx;
	assert(distmx_i != 0 && distmx_j != 0);

	const float revw = flat_params::m_rev_w;
	if (revw > 0)
		{
		uint ncol_rev, lo_i_rev, lo_j_rev;
		const uint8_t *rev_prof_i = cd_i->m_mega_prof_rev;
		asserta(rev_prof_i != 0);
		float rev_score = sw_flat_pssm(
			TD.m_scratch_rows, TD.m_TB, TD.m_scratch_pssms,
			rev_prof_i, L_i,
			pssm_j, L_j, 
			flat_alphas::m_feature_block_offsets,
			flat_alphas::m_nfeat,
			-flat_params::m_open, 
			-flat_params::m_ext,
			lo_i_rev, lo_j_rev, TD.m_path_buffer, ncol_rev);
		score -= revw*rev_score;
		}

	const float selfw = flat_params::m_self_w;
	if (selfw > 0)
		score -= selfw*(m_self_rev_scores[i] + m_self_rev_scores[j])/2;

	if (flat_params::m_lddt_w > 0)
		{
		float lddt = flat_getlddt_muscle_some_floats3(
			label_i, label_j, path,
			lo_i, L_i, lo_j, L_j, distmx_i, distmx_j);
		score += flat_params::m_lddt_w*lddt*500;
		}

	if (flat_params::m_lddtx_w > 0)
		{
		float L = (L_i + L_j)/2.0f + 50;
		float Lfactor = float(ncol)/L;

		float lddt = flat_getlddt_muscle_some_floats3(
			label_i, label_j, path,
			lo_i, L_i, lo_j, L_j, distmx_i, distmx_j);
		score += flat_params::m_lddtx_w*lddt*500*Lfactor;
		}

	if (flat_params::m_dali_w > 0)
		{
		float dali = flat_get_dali(
			label_i, label_j, path,
			lo_i, L_i, lo_j, L_j, distmx_i, distmx_j);
		score += flat_params::m_dali_w*dali*10;
		}

	if (flat_params::m_dalix_w > 0)
		{
		float dalix = flat_get_dalix(
			label_i, label_j, path,
			lo_i, L_i, lo_j, L_j,
			distmx_i, distmx_j, TD.m_colscores);
		score += flat_params::m_dalix_w*dalix*10;
		}

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
	chain_data::update_pssms(m_cdvec, m_look->get_ndom());
	if (flat_params::need_self())
		set_self_rev_scores();
	if (flat_params::need_nu_self())
		set_nu_self_rev_scores();
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
	FB.load_chains(chains);
	FB.update_params(param_names, param_values);

	flat_alphas::logme();
	flat_params::logme();

	uint nthread = GetRequestedThreadCount();
	thread_affinity ta;
	bool pin = opt(no_thread_pin) ? false : ta.shouldPin(nthread);
	FB.search(nthread, pin);
	FB.SetScoreOrder();
	FB.Bench();
	FB.WriteHits(opt(output), opt(include_self), opt(triangle));

	double align_count = double(FB.m_aln_count);
	double mu_fwd_reject_count= double(FB.m_mu_fwd_reject_count);
	double mu_combined_reject_count= double(FB.m_mu_combined_reject_count);
	ProgressLog("Mu filter fwd %.1f%%, combined %.1f%%\n",
		GetPct(mu_fwd_reject_count, align_count),
		GetPct(mu_combined_reject_count, align_count));
	}
