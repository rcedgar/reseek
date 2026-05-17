#include "myutils.h"
#include "flat_bench2.h"
#include "thread_affinity.h"
#include "flat_helpers.h"
#include "flat_alphas.h"
#include "paralign.h"
#include "cigar.h"

uint flat_bench2::m_maxL = 4000;

static FILE *s_f_nu_paths;

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

	time_t t1 = time(0);

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

	time_t t2 = time(0);
	ProgressLog("Search time %.0f secs.\n", double(t2 - t1));
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

	const uint ndom = m_look->get_ndom();
	if (m_self_rev_scores == 0)
		m_self_rev_scores = myalloc(float, ndom);

	const uint ThreadCount = GetRequestedThreadCount();

#pragma omp parallel num_threads(ThreadCount)
	{
	flat_bench2_thread_data TD(m_maxL, flat_alphas::m_nfeat);

#pragma omp for schedule(dynamic)
	for (int domidx = 0; domidx < (int)ndom; ++domidx)
		{
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
	}  // each thread destroys its TD here
	}

void flat_bench2::load_chains(const vector<flat_chain_t *> &chains)
	{
	uint nchain = uint(chains.size());
	uint ndom = m_look->get_ndom();
	asserta(nchain >= ndom);
	m_cdvec = myalloc(chain_data *, ndom);
	memset(m_cdvec, 0, ndom*sizeof(m_cdvec[0]));
	vector<flat_chain_t *> sorted_chains;
	m_look->sort_chains(chains, sorted_chains);
	chain_data::fill_chain_data_vec(sorted_chains, bits_query, m_cdvec);
	}

float flat_bench2::score_pos_pair(
	const uint8_t *mega_prof_i, uint pos_i, uint L_i,
	const uint8_t *mega_prof_j, uint pos_j, uint L_j) const
	{
	float score = 0;
	for (uint fi = 0; fi < flat_alphas::m_nfeat; ++fi)
		{
		const float *weighted_logoddsvec =
			flat_alphas::m_weighted_logoddsvec[fi];
		uint alpha_size = flat_alphas::m_alpha_sizes[fi];
		assert(weighted_logoddsvec != 0);
		uint8_t code_i = mega_prof_i[fi*L_i + pos_i];
		uint8_t code_j = mega_prof_j[fi*L_j + pos_j];
		assert(code_i < alpha_size);
		assert(code_j < alpha_size);
		score += weighted_logoddsvec[code_i*alpha_size + code_j];
		}
	return score;
	}

float flat_bench2::score_path(
	const chain_data &cd_i,
	uint lo_i,
	const chain_data &cd_j,
	uint lo_j,
	const char *path,
	uint ncol) const
	{
	const uint L_i = cd_i.m_L;
	const uint L_j = cd_j.m_L;
	uint pos_i = lo_i;
	uint pos_j = lo_j;
	bool in_gap = false;
	const float open = -flat_params::m_open;
	const float ext = -flat_params::m_ext;
	asserta(open <= 0);
	asserta(ext <= 0);
	float score= 0;
	const uint8_t *mega_prof_i = cd_i.m_mega_prof;
	const uint8_t *mega_prof_j = cd_j.m_mega_prof;
	for (uint col = 0; col < ncol; ++col)
		{
		char c = path[col];
		if (c == 'M')
			{
			assert(pos_i < L_i);
			in_gap = false;
			score += score_pos_pair(
				cd_i.m_mega_prof, pos_i, L_i,
				cd_j.m_mega_prof, pos_j, L_j);
			++pos_i;
			++pos_j;
			}
		else if (c == 'D')
			++pos_i;
		else if (c == 'I')
			++pos_j;

		if (c == 'D' || c == 'I')
			{
			if (in_gap)
				score += ext;
			else
				{
				score += open;
				in_gap = true;
				}
			}
		}
	return score;
	}

void flat_bench2::align_pair_nu_paths(
	uint pairidx, flat_bench2_thread_data &TD)
	{
	asserta(s_f_nu_paths);
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

	const int open = Paralign::m_Open;
	const int ext = Paralign::m_Ext;

	if (TD.m_parasail_result != 0)
		{
		parasail_result_free(TD.m_parasail_result);
		TD.m_parasail_result = 0;
		}

	parasail_profile_t *prof_i = cd_i->m_parasail_prof;
	asserta(prof_i != 0);

	const uint8_t *codeseq_nu_i = cd_i->m_codeseq_nu;
	const uint8_t *codeseq_nu_j = cd_j->m_codeseq_nu;
	const uint8_t *codeseq_nu_j_rev = cd_j->m_codeseq_nu_rev;
	asserta(codeseq_nu_i != 0);
	asserta(codeseq_nu_j != 0);
	asserta(codeseq_nu_j_rev != 0);

	string cigar, cigar_rev;
	uint Lo_i, Lo_j, Lo_i_rev, Lo_j_rev;
	int fwd_score, rev_score;
	string fwd_path;
	string rev_path;
	string fwd_compact_cigar;
	string rev_compact_cigar;
	{
	TD.m_parasail_result = parasail_sw_trace_striped_profile_avx2_256_16(
		prof_i, (const char *) codeseq_nu_j, L_j, open, ext);
	asserta(!(TD.m_parasail_result->flag & PARASAIL_FLAG_SATURATED));
	fwd_score = TD.m_parasail_result->score;

	parasail_cigar_t* cig = parasail_result_get_cigar_extra(
		TD.m_parasail_result,
		(const char *) codeseq_nu_i, L_i,
		(const char *) codeseq_nu_j, L_j,
		&Paralign::m_matrix, 1, 0);

	char *cig_str = parasail_cigar_decode(cig);
	cigar = string(cig_str);
	Lo_i = (uint) cig->beg_query;
	Lo_j = (uint) cig->beg_ref;
	free(cig_str);
	parasail_cigar_free(cig);
	parasail_result_free(TD.m_parasail_result);
	TD.m_parasail_result = 0;

	ExpandParaCigar_reverseDI(cigar, fwd_path);

	int check_score_fwd = Paralign::score_nu_path(
		label_i, codeseq_nu_i, Lo_i, L_i,
		label_j, codeseq_nu_j, Lo_j, L_j,
		fwd_path);

	PathToCIGAR(fwd_path.c_str(), fwd_compact_cigar);
	}

	{
	TD.m_parasail_result = parasail_sw_trace_striped_profile_avx2_256_16(
		prof_i, (const char *) codeseq_nu_j_rev, L_j, open, ext);
	asserta(!(TD.m_parasail_result->flag & PARASAIL_FLAG_SATURATED));
	rev_score = TD.m_parasail_result->score;

	parasail_cigar_t* cig_rev = parasail_result_get_cigar_extra(
		TD.m_parasail_result,
		(const char *) codeseq_nu_i, L_i,
		(const char *) codeseq_nu_j_rev, L_j,
		&Paralign::m_matrix, 1, 0);
	char *cig_str_rev = parasail_cigar_decode(cig_rev);
	cigar_rev = string(cig_str_rev);
	Lo_i_rev = (uint) cig_rev->beg_query;
	Lo_j_rev = (uint) cig_rev->beg_ref;
	free(cig_str_rev);
	parasail_cigar_free(cig_rev);
	parasail_result_free(TD.m_parasail_result);
	TD.m_parasail_result = 0;

	ExpandParaCigar_reverseDI(cigar_rev, rev_path);

	int check_score_rev = Paralign::score_nu_path(
		label_i, codeseq_nu_i, Lo_i_rev, L_i,
		label_j, codeseq_nu_j_rev, Lo_j_rev, L_j,
		rev_path);

	PathToCIGAR(rev_path.c_str(), rev_compact_cigar);
	}

// NOTE -- sometimes parasail_cigar_decode returns
// cig_rev->beg_query=0, cig_rev->beg_ref=0 and
// a CIGAR string which begins with Ds or Is
// This is a quirk more than a bug, it's a non-
// standard way to represent local alignment
	static mutex lock;
	FILE *f = s_f_nu_paths;
	lock.lock();
	fprintf(f, "%s", label_i.c_str());
	fprintf(f, "\t%s", label_j.c_str());
	fprintf(f, "\t%u", Lo_i);
	fprintf(f, "\t%u", Lo_j);
	fprintf(f, "\t%s", fwd_compact_cigar.c_str());
	fprintf(f, "\t%d", fwd_score);
	fprintf(f, "\t%u", Lo_i_rev);
	fprintf(f, "\t%u", Lo_j_rev);
	fprintf(f, "\t%s", rev_compact_cigar.c_str());
	fprintf(f, "\t%d", rev_score);
	fprintf(f, "\n");
	lock.unlock();
	}

void flat_bench2::align_pair(
	uint pairidx, flat_bench2_thread_data &TD)
	{
	if (m_nu_paths)
		{
		align_pair_nu_paths(pairidx, TD);
		return;
		}
	m_Scores[pairidx] = 0;

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

	float nu_rev_score = 0;
	if (flat_params::m_nu_filter_min_fwd_score > 0)
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
		if (fwd_score < flat_params::m_nu_filter_min_fwd_score)
			{
			++m_mu_fwd_reject_count;
			return;
			}

		parasail_result_free(TD.m_parasail_result);
		parasail_profile_t *prof_i_rev = cd_i->m_parasail_prof_rev;
		asserta(prof_i_rev != 0);
		TD.m_parasail_result = parasail_sw_striped_profile_avx2_256_16(
			prof_i_rev, (const char *) codeseq_nu_j, L_j, open, ext);
		asserta(!(TD.m_parasail_result->flag & PARASAIL_FLAG_SATURATED));
		nu_rev_score = (float) TD.m_parasail_result->score;

		float self_score = (m_nu_self_rev_scores[i] + 
			m_nu_self_rev_scores[j])/2.0f;

		const float revw = flat_params::m_nu_filter_rev_w;
		const float selfw = flat_params::m_nu_filter_self_w;

		float nu_combined_score =
			float(fwd_score) -
			selfw*self_score -
			revw*nu_rev_score;
		if (nu_combined_score < flat_params::m_nu_filter_min_combined_score)
			{
			++m_mu_combined_reject_count;
			return;
			}
		}
	
	const uint8_t *prof_i = cd_i->m_mega_prof;
	const float *pssm_j = cd_j->m_mega_pssm;

	uint lo_i, lo_j, ncol;
	float score = 0;
	float mega_fwd_score = sw_flat_pssm(
		TD.m_scratch_rows, TD.m_TB, TD.m_scratch_pssms,
		prof_i, L_i,
		pssm_j, L_j, 
		flat_alphas::m_feature_block_offsets,
		flat_alphas::m_nfeat,
		-flat_params::m_open, 
		-flat_params::m_ext,
		lo_i, lo_j, TD.m_path_buffer, ncol);
	const string path = string(TD.m_path_buffer);
	++m_mega_fwd_test_count;
	if (mega_fwd_score < flat_params::m_mega_filter_min_fwd)
		return;
	score += mega_fwd_score;
	++m_mega_fwd_pass_count;

	//float score2 = score_path(
	//	*cd_i, lo_i, *cd_j, lo_j, TD.m_path_buffer, ncol);
	//asserta(score2 == mega_fwd_score);

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
	score += flat_params::m_nurev_w*nu_rev_score;

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
	asserta(m_cdvec != 0);
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
	const bool nu_paths = optset_output2;
	if (nu_paths)
		s_f_nu_paths = CreateStdioFile(opt(output2));

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

	string varstr2;
	flat_make_varstr(varstr2);
	Log("varstr2=\n");
	Log("%s\n", varstr2.c_str());

	vector<string> peaker_spec_lines;
	flat_make_peaker_spec(peaker_spec_lines);
	for (uint i = 0; i < uint(peaker_spec_lines.size()); ++i)
		Log("%s\n", peaker_spec_lines[i].c_str());

	uint nthread = GetRequestedThreadCount();
	thread_affinity ta;
	bool pin = opt(no_thread_pin) ? false : ta.shouldPin(nthread);
	FB.m_nu_paths = nu_paths;
	FB.search(nthread, pin);
	if (nu_paths)
		{
		CloseStdioFile(s_f_nu_paths);
		return;
		}

	FB.SetScoreOrder();
	FB.Bench();
	FB.WriteHits(opt(output), opt(include_self), opt(triangle));

	double align_count = double(FB.m_aln_count);
	double mega_fwd_test_count = double(FB.m_mega_fwd_test_count);
	double mega_fwd_pass_count = double(FB.m_mega_fwd_pass_count);
	double mu_fwd_reject_count = double(FB.m_mu_fwd_reject_count);
	double mu_combined_reject_count = double(FB.m_mu_combined_reject_count);
	double mega_passed_pct = GetPct(mega_fwd_pass_count, mega_fwd_test_count);
	ProgressLog("Mu filter fwd %.1f%%, combined %.1f%%, total %.1f%%\n",
		GetPct(mu_fwd_reject_count, align_count),
		GetPct(mu_combined_reject_count, align_count),
		GetPct(mu_fwd_reject_count+mu_combined_reject_count, align_count));
	ProgressLog("Mega fwd filter passed %.1f%%\n", mega_passed_pct);
	}
