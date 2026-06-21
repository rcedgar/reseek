#include "myutils.h"
#include "flat_bench2.h"
#include "thread_affinity.h"
#include "flat_helpers.h"
#include "flat_params.h"
#include "paralign.h"
#include "parasail_nomalloc.h"
#include "cigar.h"
#include "seqdb.h"
#include "getticks.h"

#define WRITE_NU_SELF_REV_SCORES	0
#define WRITE_TS_TERMS				1

#if WRITE_TS_TERMS
static FILE *s_fts;
static mutex s_ts_lock;
#endif


static FILE *s_f_nu_paths;
static FILE *s_f_mega_paths;

void flat_bench2::search(uint nthread, bool pin_threads)
	{
	FastBench::Alloc();

	const uint NQ = SIZE(m_Labels);
	const uint PairCount = triangle_get_K(NQ);

	ProgressStep(0, PairCount, "Aligning");
	m_next_pairidx = 0;
	m_aln_count = 0;
	m_nu_fwd_reject_count = 0;
	m_nu_combined_reject_count = 0;
	m_nu_pass_count = 0;

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

	time_t t2 = time(0);
	ProgressStep(PairCount-1, PairCount,
		"Search time %.0f secs", double(t2 - t1));
	Log("Search time %.0f secs\n", double(t2 - t1));
	}

void flat_bench2::static_thread_body(flat_bench2 *FB, uint threadidx)
	{
	FB->thread_body(threadidx);
	}

void flat_bench2::static_thread_body_set_mega_self_rev_scores(flat_bench2 *FB, uint threadidx)
	{
	FB->thread_body_set_mega_self_rev_scores(threadidx);
	}

void flat_bench2::thread_body_set_mega_self_rev_scores(uint threadidx)
	{
	const uint ndom = m_look->get_ndom();

	uint nfeat = m_params->m_nfeat;
	asserta(nfeat > 0);

	flat_bench2_thread_data TD(nfeat);
	for (;;)
		{
		uint domidx = m_next_domidx++;
		if (domidx >= ndom)
			return;
		set_mega_self_rev_score(domidx, TD);
		}
	}

void flat_bench2::thread_body(uint threadidx)
	{
	const uint NQ = SIZE(m_Labels);
	const uint npair = triangle_get_K(NQ);

	uint nfeat = m_params->m_nfeat;
	asserta(nfeat > 0);

	flat_bench2_thread_data TD(nfeat);
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
	if (!m_params->need_nu_self())
		{
		asserta(m_nu_self_rev_scores == 0);
		return;
		}

	flat_bench2_thread_data TD(m_params->m_nfeat);
	const uint ndom = m_look->get_ndom();
	if (m_nu_self_rev_scores == 0)
		m_nu_self_rev_scores = myalloc(float, ndom);

#if WRITE_NU_SELF_REV_SCORES
	FILE *ftmp = CreateStdioFile("nu_self_rev_scores.tmp");
#endif
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		const chain_data *cd = m_cdvec[domidx];
		const int open = Paralign::m_Open;
		const int ext = Paralign::m_Ext;

		if (TD.m_parasail_result != 0)
			parasail_result_free(TD.m_parasail_result);
		parasail_profile_t *prof_rev = cd->m_parasail_prof_rev;
		asserta(prof_rev != 0);
		const uint8_t *codeseq_nu = cd->m_codeseq_nu;
		asserta(codeseq_nu != 0);
		int self_rev_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
			prof_rev, (const char *) codeseq_nu, cd->m_L, open, ext,
			TD.m_parasail_nomalloc_workspace, TD.m_parasail_nomalloc_workspace_bytes,
			0, 0, 0);
		m_nu_self_rev_scores[domidx] = float(self_rev_score);
#if WRITE_NU_SELF_REV_SCORES
		fprintf(ftmp, "%s\t%d\n", m_look->get_dom(domidx).c_str(), self_rev_score);
#endif

		}
#if WRITE_NU_SELF_REV_SCORES
	CloseStdioFile(ftmp);
#endif
	}

void flat_bench2::set_mega_self_rev_score(
	uint domidx, flat_bench2_thread_data &TD)
	{
	const chain_data *cd = m_cdvec[domidx];
	const uint8_t *prof = cd->m_mega_prof;
	const float *pssm = cd->m_mega_pssm_rev;
	asserta(prof != 0);
	asserta(pssm != 0);
	const uint L = cd->m_L;
	uint lo_i, lo_j, ncol;
	float score = sw_flat_pssm(
		TD.m_scratch_rows, TD.m_TB, TD.m_scratch_pssms,
		prof, L, pssm, L,
		m_params->m_feature_block_offsets,
		m_params->m_nfeat,
		-m_params->m_open,
		-m_params->m_ext,
		lo_i, lo_j, TD.m_path_buffer, ncol);

	asserta(!isinf(score));
	asserta(!isnan(score));
	m_self_rev_scores[domidx] = score;
	}

void flat_bench2::set_mega_self_rev_scores()
	{
	const uint nthread = GetRequestedThreadCount();
	const uint ndom = m_look->get_ndom();
	if (m_self_rev_scores == 0)
		m_self_rev_scores = myalloc(float, ndom);
	memset(m_self_rev_scores, 0xff, ndom*sizeof(float));

	m_next_domidx = 0;
	vector<thread *> ts;
	for (uint threadidx = 0; threadidx < nthread; ++threadidx)
		{
		thread *t = new thread(static_thread_body_set_mega_self_rev_scores, this, threadidx);
		ts.push_back(t);
		}
	for (uint threadidx = 0; threadidx < nthread; ++threadidx)
		ts[threadidx]->join();
	for (uint threadidx = 0; threadidx < nthread; ++threadidx)
		delete ts[threadidx];

	for (uint i = 0; i < ndom; ++i)
		{
		float score = m_self_rev_scores[i];
		asserta(!isnan(score));
		asserta(!isinf(score));
		}
	}

void flat_bench2::write_nu_hexfasta(const string &fn)
	{
	if (fn == "") return;
	const uint ndom = m_look->get_ndom();
	FILE *f = CreateStdioFile(fn);
	for (uint i = 0; i < ndom; ++i)
		{
		const chain_data *cd = m_cdvec[i];
		const string &label = cd->m_label;
		const uint8_t *codeseq_nu = cd->m_codeseq_nu;
		codeseq_to_hexfasta(f, label, codeseq_nu, cd->m_L);
		}
	CloseStdioFile(f);
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
	chain_data::fill_chain_data_vec(
		*m_params, sorted_chains, bits_query, m_cdvec);
	}

float flat_bench2::score_pos_pair(
	const uint8_t *mega_prof_i, uint pos_i, uint L_i,
	const uint8_t *mega_prof_j, uint pos_j, uint L_j) const
	{
	float score = 0;
	for (uint fi = 0; fi < m_params->m_nfeat; ++fi)
		{
		const float *weighted_logoddsvec =
			m_params->m_weighted_logoddsvec[fi];
		uint alpha_size = m_params->m_alpha_sizes[fi];
		assert(weighted_logoddsvec != 0);
		uint8_t code_i = mega_prof_i[fi*L_i + pos_i];
		uint8_t code_j = mega_prof_j[fi*L_j + pos_j];
		assert(code_i < alpha_size);
		assert(code_j < alpha_size);
		score += weighted_logoddsvec[code_i*alpha_size + code_j];
		}
	return score;
	}

float flat_bench2::calc_ts(
	uint i, uint j,
	const chain_data &cd_i,
	const chain_data &cd_j,
	uint fwd_lo_i, uint fwd_lo_j,
	const char *fwd_path,
	uint fwd_ncol,
	uint rev_lo_i, uint rev_lo_j,
	const char *rev_path,
	uint rev_ncol,
	flat_bench2_thread_data &TD)
	{
	//const string &label_i = cd_i.m_label;
	//const string &label_j = cd_j.m_label;

	const uint L_i = cd_i.m_L;
	const uint L_j = cd_j.m_L;

	const sid_t *distmx_i = cd_i.m_distmx;
	const sid_t *distmx_j = cd_j.m_distmx;

	float ts = 0;
	const float revw = m_params->m_rev_w;
	if (revw > 0)
		{
		float score_rev = score_path(
			cd_i, rev_lo_i,
			cd_j, rev_lo_j,
			rev_path, rev_ncol);
		ts -= revw*score_rev;
		}

	uint nmatch = path2posvecs3(TD.m_path_buffer, fwd_ncol,
		fwd_lo_i, L_i, fwd_lo_j, L_j, TD.m_pos_is, TD.m_pos_js, flat_params::m_maxL);

	float nu_rev_score = 0;//@@TODO
	asserta(m_params->m_nurev_w == 0);
	ts += m_params->m_nurev_w*nu_rev_score;

	const float selfw = m_params->m_self_w;
	if (selfw > 0)
		ts -= selfw*(m_self_rev_scores[i] + m_self_rev_scores[j])/2;

	if (m_params->m_lddt_w > 0)
		{
		float lddt = flat_getlddt_muscle_some_floats(
			TD.m_pos_is, L_i,
			TD.m_pos_js, L_j,
			nmatch, distmx_i, distmx_j,
			TD.m_considered_vec, TD.m_preserved_vec);
		ts += m_params->m_lddt_w*lddt*500;
		}

	if (m_params->m_lddtx_w > 0)
		{
		asserta(false);//lddtx screwed up?
		float L = (L_i + L_j)/2.0f + 50;
		float Lfactor = float(fwd_ncol)/L;

		float lddt = flat_getlddt_muscle_some_floats( // TODO this is not lddtx!?
			TD.m_pos_is, L_i,
			TD.m_pos_js, L_j,
			nmatch, distmx_i, distmx_j,
			TD.m_considered_vec, TD.m_preserved_vec);
		ts += m_params->m_lddtx_w*lddt*500*Lfactor;
		}

	if (m_params->m_dali_w > 0)
		{
		float dali = flat_get_dali4(
			TD.m_pos_is, L_i,
			TD.m_pos_js, L_j,
			nmatch, distmx_i, distmx_j);
		ts += m_params->m_dali_w*dali*10;
		}

	if (m_params->m_dalix_w > 0)
		{
		Die("TODO");
		}

	return ts;
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
	const float open = -m_params->m_open;
	const float ext = -m_params->m_ext;
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

void flat_bench2::align_pair_input_mega_paths(
	uint pairidx, flat_bench2_thread_data &TD)
	{
	const string &fwd_cigar = m_mega_path_cigar_fwds[pairidx];
	if (fwd_cigar.empty())
		{
		m_Scores[pairidx] = 0;
		return;
		}

	uint i, j;
	uint NQ = uint(m_Labels.size());
	triangle_k_to_ij(pairidx, NQ, i, j);

	const string &rev_cigar = m_mega_path_cigar_revs[pairidx];

	uint fwd_lo_i = m_mega_path_lo_i_fwds[pairidx];
	uint fwd_lo_j = m_mega_path_lo_j_fwds[pairidx];

	uint rev_lo_i = m_mega_path_lo_i_revs[pairidx];
	uint rev_lo_j = m_mega_path_lo_j_revs[pairidx];

	string fwd_path, rev_path;
	CIGARToPath(fwd_cigar, fwd_path, true);
	CIGARToPath(rev_cigar, rev_path, true);

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
	asserta(L_i <= flat_params::m_maxL); // TODO=maxL
	asserta(L_j <= flat_params::m_maxL); // TODO=maxL

	//float mega_score = score_path(
	//	*cd_i, lo_i,
	//	*cd_j, lo_j,
	//	path.c_str(), uint(path.size()));

	float ts = calc_ts(
		i, j, *cd_i, *cd_j,
		fwd_lo_i, fwd_lo_j,
		fwd_path.c_str(), uint(fwd_path.size()),
		rev_lo_i, rev_lo_j,
		rev_path.c_str(), uint(rev_path.size()),
		TD);

	asserta(!isnan(ts));
	asserta(!isinf(ts));
	m_Scores[pairidx] = ts;
	}

void flat_bench2::align_pair_output_nu_paths(
	uint pairidx, flat_bench2_thread_data &TD)
	{
	const float MIN_FWD_SCORE = float(optset_minscore ? opt(minscore) : 120);

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
	asserta(L_i <= flat_params::m_maxL); // TODO=maxL
	asserta(L_j <= flat_params::m_maxL); // TODO=maxL

	const int open = Paralign::m_Open;
	const int ext = Paralign::m_Ext;

	if (TD.m_parasail_result != 0)
		{
		parasail_result_free(TD.m_parasail_result);
		TD.m_parasail_result = 0;
		}

	parasail_profile_t *prof_i = cd_i->m_parasail_prof;
	parasail_profile_t *prof_i_rev = cd_i->m_parasail_prof_rev;
	asserta(prof_i != 0);
	asserta(prof_i_rev != 0);

	const uint8_t *codeseq_nu_i = cd_i->m_codeseq_nu;
	const uint8_t *codeseq_nu_j = cd_j->m_codeseq_nu;
	asserta(codeseq_nu_i != 0);
	asserta(codeseq_nu_j != 0);

	string para_cigar_fwd;
	string para_cigar_rev;
	uint lo_i_fwd = UINT_MAX;
	uint lo_j_fwd = UINT_MAX;
	uint lo_i_rev = UINT_MAX;
	uint lo_j_rev = UINT_MAX;
	int score_fwd = 0;
	int score_rev = 0;
	string path_fwd;
	string path_rev;
	string compact_cigar_fwd;
	string compact_cigar_rev;
	{
	if (TD.m_parasail_result != 0)
		parasail_result_free(TD.m_parasail_result);
	TD.m_parasail_result = parasail_sw_trace_striped_profile_avx2_256_16(
		prof_i, (const char *) codeseq_nu_j, L_j, open, ext);
	asserta(!(TD.m_parasail_result->flag & PARASAIL_FLAG_SATURATED));
	score_fwd = TD.m_parasail_result->score;
	if (score_fwd < MIN_FWD_SCORE)
		return;

	parasail_cigar_t* cig = parasail_result_get_cigar_extra(
		TD.m_parasail_result,
		(const char *) codeseq_nu_i, L_i,
		(const char *) codeseq_nu_j, L_j,
		&Paralign::m_matrix, 1, 0);

	char *cig_str = parasail_cigar_decode(cig);
	para_cigar_fwd = string(cig_str);
	uint para_lo_i_fwd = (uint) cig->beg_query;
	uint para_lo_j_fwd = (uint) cig->beg_ref;
	free(cig_str);
	parasail_cigar_free(cig);

	string path_fwd2;
	uint lo_i_fwd2, lo_j_fwd2;
	parasail_result_to_path(TD.m_parasail_result,
		L_i, L_j, lo_i_fwd2, lo_j_fwd2, path_fwd2);

	parasail_result_free(TD.m_parasail_result);
	TD.m_parasail_result = 0;

	parasail_cigar_to_path(para_cigar_fwd,
		para_lo_i_fwd, para_lo_j_fwd,
		lo_i_fwd, lo_j_fwd, path_fwd);

	asserta(lo_i_fwd2 == lo_i_fwd);
	asserta(lo_j_fwd2 == lo_j_fwd);
	asserta(path_fwd2 == path_fwd);

	int check_score_fwd = Paralign::score_nu_path(
		label_i, codeseq_nu_i, lo_i_fwd, L_i,
		label_j, codeseq_nu_j, lo_j_fwd, L_j,
		path_fwd);
	if (check_score_fwd != score_fwd)
		{
		Log("i=%u=%s\n", i, label_i.c_str());
		Log("j=%u=%s\n", j, label_i.c_str());
		Log("para_lo_i_fwd=%u\n", para_lo_i_fwd);
		Log("para_lo_j_fwd=%u\n", para_lo_j_fwd);
		Log("para_cigar_fwd=%s\n", para_cigar_fwd.c_str());
		Log("path_fwd=%s\n", path_fwd.c_str());
		Warning("check_score_fwd=%d != score_fwd=%d",
			check_score_fwd, score_fwd);
		}

	PathToCIGAR(path_fwd.c_str(), compact_cigar_fwd);
	}

	{
	TD.m_parasail_result = parasail_sw_trace_striped_profile_avx2_256_16(
		prof_i_rev, (const char *) codeseq_nu_j, L_j, open, ext);
	asserta(!(TD.m_parasail_result->flag & PARASAIL_FLAG_SATURATED));
	score_rev = TD.m_parasail_result->score;

	//parasail_cigar_t* cig_rev = parasail_result_get_cigar_extra(
	//	TD.m_parasail_result,
	//	(const char *) codeseq_nu_i, L_i,
	//	(const char *) codeseq_nu_j_rev, L_j,
	//	&Paralign::m_matrix, 1, 0);
	//char *cig_str_rev = parasail_cigar_decode(cig_rev);
	//para_cigar_rev = string(cig_str_rev);
	//uint para_lo_i_rev = (uint) cig_rev->beg_query;
	//uint para_lo_j_rev = (uint) cig_rev->beg_ref;
	//free(cig_str_rev);
	//parasail_cigar_free(cig_rev);

	string path_rev2;
	uint lo_i_rev2, lo_j_rev2;
	parasail_result_to_path(TD.m_parasail_result,
		L_i, L_j, lo_i_rev2, lo_j_rev2, path_rev2);

	parasail_result_free(TD.m_parasail_result);
	TD.m_parasail_result = 0;

	//parasail_cigar_to_path(para_cigar_rev,
	//	para_lo_i_rev, para_lo_j_rev,
	//	lo_i_rev, lo_j_rev, path_rev);

	asserta(lo_i_rev2 == lo_i_rev);
	asserta(lo_j_rev2 == lo_j_rev);
	asserta(path_rev2 == path_rev);

	//int check_score_rev = Paralign::score_nu_path(
	//	label_i, codeseq_nu_i, lo_i_rev, L_i,
	//	label_j, codeseq_nu_j_rev, lo_j_rev, L_j,
	//	path_rev);

	PathToCIGAR(path_rev.c_str(), compact_cigar_rev);
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
	fprintf(f, "\t%u", lo_i_fwd);
	fprintf(f, "\t%u", lo_j_fwd);
	fprintf(f, "\t%s", compact_cigar_fwd.c_str());
	fprintf(f, "\t%d", score_fwd);
	fprintf(f, "\t%u", lo_i_rev);
	fprintf(f, "\t%u", lo_j_rev);
	fprintf(f, "\t%s", compact_cigar_rev.c_str());
	fprintf(f, "\t%d", score_rev);
	fprintf(f, "\n");
	lock.unlock();
	}

void flat_bench2::align_pair_nu_only(
	uint pairidx, flat_bench2_thread_data &TD)
	{
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
	asserta(L_i <= flat_params::m_maxL); // TODO=maxL
	asserta(L_j <= flat_params::m_maxL); // TODO=maxL

	int nu_rev_score = 0;
	const int open = Paralign::m_Open;
	const int ext = Paralign::m_Ext;

	if (TD.m_parasail_result != 0)
		parasail_result_free(TD.m_parasail_result);
	parasail_profile_t *prof_i = cd_i->m_parasail_prof;
	asserta(prof_i != 0);
	const uint8_t *codeseq_nu_j = cd_j->m_codeseq_nu;
	int fwd_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
		prof_i, (const char *) codeseq_nu_j, L_j, open, ext,
		TD.m_parasail_nomalloc_workspace, TD.m_parasail_nomalloc_workspace_bytes);
	if (fwd_score < m_params->m_nu_filter_min_fwd_score)
		{
		++m_nu_fwd_reject_count;
		return;
		}

	parasail_profile_t *prof_i_rev = cd_i->m_parasail_prof_rev;
	asserta(prof_i_rev != 0);
	nu_rev_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
		prof_i_rev, (const char *) codeseq_nu_j, L_j, open, ext,
		TD.m_parasail_nomalloc_workspace, TD.m_parasail_nomalloc_workspace_bytes);

	float self_score = (m_nu_self_rev_scores[i] + 
		m_nu_self_rev_scores[j])/2.0f;

	const float revw = m_params->m_nu_filter_rev_w;
	const float selfw = m_params->m_nu_filter_self_w;

	float nu_combined_score =
		float(fwd_score) -
		selfw*self_score -
		revw*nu_rev_score;
	if (nu_combined_score < m_params->m_nu_filter_min_combined_score)
		{
		++m_nu_combined_reject_count;
		return;
		}
	++m_nu_pass_count;
	asserta(!isnan(nu_combined_score));
	asserta(!isinf(nu_combined_score));
	m_Scores[pairidx] = nu_combined_score;
	}

void flat_bench2::align_pair_timealn(
	uint pairidx, flat_bench2_thread_data &TD)
	{
	static mutex lock;
	lock.lock();
	static uint counter = 0;
	static TICKS sumticks_nu_scoreonly = 0;
	static TICKS sumticks_nu_path = 0;
	static TICKS sumticks_mega_scoreonly = 0;
	static TICKS sumticks_mega_path = 0;

	if (counter > 0 && counter%10000 == 0)
		{
		double nu_scoreonly = double(sumticks_nu_scoreonly);
		double nu_path = double(sumticks_nu_path);
		double mega_scoreonly = double(sumticks_mega_scoreonly);
		double mega_path = double(sumticks_mega_path);

		ProgressLog("\n");
		ProgressLog("nu_scoreonly    %8.3g  %6.2f x\n", nu_scoreonly, mega_path/nu_scoreonly);
		ProgressLog("nu_path         %8.3g  %6.2f x\n", nu_path, mega_path/nu_path);
		ProgressLog("mega_scoreonly  %8.3g  %6.2f x\n", mega_scoreonly, mega_path/mega_scoreonly);
		ProgressLog("mega_path       %8.3g  %6.2f x\n", mega_path, 1.0);
		ProgressLog("\n");
		if (counter == 40000)
			Die("timealn done");
		}

	++counter;
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
	asserta(L_i <= flat_params::m_maxL); // TODO=maxL
	asserta(L_j <= flat_params::m_maxL); // TODO=maxL

	const int open = Paralign::m_Open;
	const int ext = Paralign::m_Ext;

	TICKS t1 = GetClockTicks();
	if (TD.m_parasail_result != 0)
		{
		parasail_result_free(TD.m_parasail_result);
		TD.m_parasail_result = 0;
		}
	parasail_profile_t *para_prof_i = cd_i->m_parasail_prof;
	asserta(para_prof_i != 0);
	const uint8_t *codeseq_nu_j = cd_j->m_codeseq_nu;
	parasail_sw_striped_profile_avx2_256_16_nomalloc(
		para_prof_i, (const char *) codeseq_nu_j, L_j, open, ext,
		TD.m_parasail_nomalloc_workspace, TD.m_parasail_nomalloc_workspace_bytes);

	TICKS t2 = GetClockTicks();
	sumticks_nu_scoreonly += (t2 - t1);

	if (TD.m_parasail_result != 0)
		{
		parasail_result_free(TD.m_parasail_result);
		TD.m_parasail_result = 0;
		}
	TD.m_parasail_result = parasail_sw_trace_striped_profile_avx2_256_16(
		para_prof_i, (const char *) codeseq_nu_j, L_j, open, ext);

	TICKS t3 = GetClockTicks();
	sumticks_nu_path += (t3 - t2);

	const uint8_t *mega_prof_i = cd_i->m_mega_prof;
	const float *pssm_j = cd_j->m_mega_pssm;

	uint lo_i, lo_j, ncol;
	sw_flat_pssm(
		TD.m_scratch_rows, TD.m_TB, TD.m_scratch_pssms,
		mega_prof_i, L_i,
		pssm_j, L_j, 
		m_params->m_feature_block_offsets,
		m_params->m_nfeat,
		-m_params->m_open, 
		-m_params->m_ext,
		lo_i, lo_j, TD.m_path_buffer, ncol);
	TICKS t4 = GetClockTicks();
	sumticks_mega_path += (t4 - t3);

	sw_flat_pssm_scoreonly(
		TD.m_scratch_rows, TD.m_scratch_pssms,
		mega_prof_i, L_i,
		pssm_j, L_j, 
		m_params->m_feature_block_offsets,
		m_params->m_nfeat,
		-m_params->m_open, 
		-m_params->m_ext);
	TICKS t5 = GetClockTicks();
	sumticks_mega_scoreonly += (t5 - t4);

	lock.unlock();
	}

void flat_bench2::align_pair_single_feature(uint pairidx)
	{
	float sw_flatmx_scoreonly(XDPMem &Mem,
		const uint8_t *A, uint LA,
		const uint8_t *B, uint LB,
		const float *flatmx, uint alpha_size,
		float Open, float Ext);

	asserta(m_feature_logodds != 0);

	uint NQ = uint(m_Labels.size());
	uint i, j;
	triangle_k_to_ij(pairidx, NQ, i, j);
	const uint8_t *codeseq_i = m_feature_codeseq_vec[i];
	const uint8_t *codeseq_j = m_feature_codeseq_vec[j];
	uint L_i = m_feature_codeseq_lengths[i];
	uint L_j = m_feature_codeseq_lengths[j];

	thread_local XDPMem Mem;
	float score = sw_flatmx_scoreonly(
		Mem,
		codeseq_i, L_i,
		codeseq_j, L_j,
		m_feature_logodds, m_feature_alpha_size,
		m_feature_gap_open, m_feature_gap_ext);
	m_Scores[pairidx] = score;
	}

void flat_bench2::align_pair(
	uint pairidx, flat_bench2_thread_data &TD)
	{
	if (m_single_feature)
		{
		align_pair_single_feature(pairidx);
		return;
		}
	if (m_timealn)
		{
		align_pair_timealn(pairidx, TD);
		return;
		}
	if (m_nu_only)
		{
		align_pair_nu_only(pairidx, TD);
		return;
		}
	if (m_output_nu_paths)
		{
		align_pair_output_nu_paths(pairidx, TD);
		return;
		}
	else if (m_input_mega_paths)
		{
		align_pair_input_mega_paths(pairidx, TD);
		return;
		}

	m_Scores[pairidx] = 0;

	uint NQ = uint(m_Labels.size());
	uint i, j;
	triangle_k_to_ij(pairidx, NQ, i, j);

	const chain_data *cd_i = m_cdvec[i];
	const chain_data *cd_j = m_cdvec[j];

	const uint L_i = cd_i->m_L;
	const uint L_j = cd_j->m_L;
	asserta(L_i <= flat_params::m_maxL); // TODO=maxL
	asserta(L_j <= flat_params::m_maxL); // TODO=maxL

	float nu_rev_score = 0;
	int nu_fwd_hi_i = -1;
	int nu_fwd_hi_j = -1;
	if (m_params->m_nu_filter_min_fwd_score > 0)
		{
		const int open = Paralign::m_Open;
		const int ext = Paralign::m_Ext;

		parasail_profile_t *prof_i = cd_i->m_parasail_prof;
		asserta(prof_i != 0);
		const uint8_t *codeseq_nu_j = cd_j->m_codeseq_nu;
		int fwd_score = parasail_sw_striped_profile_avx2_256_16_nomalloc(
			prof_i, (const char *) codeseq_nu_j, L_j, open, ext,
			TD.m_parasail_nomalloc_workspace, TD.m_parasail_nomalloc_workspace_bytes,
			&nu_fwd_hi_i, &nu_fwd_hi_j);
		if (fwd_score < m_params->m_nu_filter_min_fwd_score)
			{
			++m_nu_fwd_reject_count;
			return;
			}

		parasail_profile_t *prof_i_rev = cd_i->m_parasail_prof_rev;
		asserta(prof_i_rev != 0);
		nu_rev_score = (float) parasail_sw_striped_profile_avx2_256_16_nomalloc(
			prof_i_rev, (const char *) codeseq_nu_j, L_j, open, ext,
			TD.m_parasail_nomalloc_workspace, TD.m_parasail_nomalloc_workspace_bytes);

		float self_score = (m_nu_self_rev_scores[i] + 
			m_nu_self_rev_scores[j])/2.0f;

		const float revw = m_params->m_nu_filter_rev_w;
		const float selfw = m_params->m_nu_filter_self_w;

		float nu_combined_score =
			float(fwd_score) -
			selfw*self_score -
			revw*nu_rev_score;
		if (nu_combined_score < m_params->m_nu_filter_min_combined_score)
			{
			++m_nu_combined_reject_count;
			return;
			}
		++m_nu_pass_count;
		}
	
	const uint8_t *prof_i = cd_i->m_mega_prof;
	const float *pssm_j = cd_j->m_mega_pssm;
	const float *pssm_j_rev = cd_j->m_mega_pssm_rev;

	uint lo_i, lo_j, ncol;
	float score = 0;
	float mega_fwd_score = sw_flat_pssm(
		TD.m_scratch_rows, TD.m_TB, TD.m_scratch_pssms,
		prof_i, L_i,
		pssm_j, L_j, 
		m_params->m_feature_block_offsets,
		m_params->m_nfeat,
		-m_params->m_open, 
		-m_params->m_ext,
		lo_i, lo_j, TD.m_path_buffer, ncol);
	//const string path = string(TD.m_path_buffer);
	++m_mega_fwd_test_count;
	if (mega_fwd_score < m_params->m_mega_filter_min_fwd)
		return;
	asserta(!isnan(mega_fwd_score));
	asserta(!isinf(mega_fwd_score));
	asserta(!isnan(score));
	asserta(!isinf(score));
	score += mega_fwd_score;
	++m_mega_fwd_pass_count;

#if WRITE_TS_TERMS
	string ts_str;
	Psa(ts_str, "%s\t%s", cd_i->m_label.c_str(), cd_j->m_label.c_str());
	Psa(ts_str, "\t%.3g", mega_fwd_score);
#endif

	uint nmatch = path2posvecs3(TD.m_path_buffer, ncol,
		lo_i, L_i, lo_j, L_j, TD.m_pos_is, TD.m_pos_js, flat_params::m_maxL);

	if (s_f_mega_paths != 0)
		{
		float score2 = score_path(
			*cd_i, lo_i, *cd_j, lo_j, TD.m_path_buffer, ncol);
		asserta(score2 == mega_fwd_score);
		string cigar;
		PathToCIGAR(TD.m_path_buffer, cigar);

		uint lo_i_rev, lo_j_rev, ncol_rev;
		float mega_score_rev = sw_flat_pssm(
			TD.m_scratch_rows, TD.m_TB, TD.m_scratch_pssms,
			prof_i, L_i,
			pssm_j_rev, L_j, 
			m_params->m_feature_block_offsets,
			m_params->m_nfeat,
			-m_params->m_open, 
			-m_params->m_ext,
			lo_i_rev, lo_j_rev, TD.m_path_buffer, ncol_rev);
		string cigar_rev;
		PathToCIGAR(TD.m_path_buffer, cigar_rev);

		uint sfidx_i = m_look->m_domidx2sfidx[i];
		uint sfidx_j = m_look->m_domidx2sfidx[j];
		string &sf_i = m_look->m_sfs[sfidx_i];
		string &sf_j = m_look->m_sfs[sfidx_j];

		string label_i = cd_i->m_label;
		string label_j = cd_j->m_label;
		trunc_label(label_i);
		trunc_label(label_j);
		asserta(label_i == m_look->get_dom(i));
		asserta(label_j == m_look->get_dom(j));

		uint cl_i = UINT_MAX;
		uint cl_j = UINT_MAX;
		uint d = find_closest_point(
			cigar, lo_i, lo_j, L_i, L_j,
			nu_fwd_hi_i, nu_fwd_hi_j,
			cl_i, cl_j);
	
		static mutex lock;
		FILE *f = s_f_mega_paths;
		lock.lock();
		fprintf(f, "%s/%s", label_i.c_str(), sf_i.c_str());  // 0
		fprintf(f, "\t%s/%s", label_j.c_str(), sf_j.c_str());  // 1
		fprintf(f, "\t%u", L_i);  // 2
		fprintf(f, "\t%u", L_j);  // 3
		fprintf(f, "\t%u", lo_i);  // 4
		fprintf(f, "\t%u", lo_j);  // 5
		fprintf(f, "\t%u", lo_i_rev);  // 6
		fprintf(f, "\t%u", lo_j_rev);  // 7
		fprintf(f, "\t%s", cigar.c_str());  // 8
		fprintf(f, "\t%s", cigar_rev.c_str());  // 9
		fprintf(f, "\t%d", nu_fwd_hi_i);  // 10
		fprintf(f, "\t%d", nu_fwd_hi_j);  // 11
		fprintf(f, "\t%u", cl_i);  // 12
		fprintf(f, "\t%u", cl_j);  // 13
		fprintf(f, "\t%u", d);  // 14
		fprintf(f, "\t%.1f", mega_fwd_score);  // 15
		fprintf(f, "\t%.1f", mega_score_rev);  // 16
		fprintf(f, "\n");
		lock.unlock();
		}

	const sid_t *distmx_i = cd_i->m_distmx;
	const sid_t *distmx_j = cd_j->m_distmx;
	assert(distmx_i != 0 && distmx_j != 0);

	const float revw = m_params->m_rev_w;
	if (revw > 0)
		{
		uint ncol_rev, lo_i_rev, lo_j_rev;
		float score_rev = sw_flat_pssm(
			TD.m_scratch_rows, TD.m_TB, TD.m_scratch_pssms,
			prof_i, L_i,
			pssm_j_rev, L_j, 
			m_params->m_feature_block_offsets,
			m_params->m_nfeat,
			-m_params->m_open, 
			-m_params->m_ext,
			lo_i_rev, lo_j_rev, TD.m_path_buffer, ncol_rev);
		score -= revw*score_rev;
		asserta(!isnan(score_rev));
		asserta(!isinf(score_rev));
		asserta(!isnan(score));
		asserta(!isinf(score));
#if WRITE_TS_TERMS
		Psa(ts_str, "\t%.3g", score_rev);
#endif
		}
	score += m_params->m_nurev_w*nu_rev_score;
#if WRITE_TS_TERMS
	Psa(ts_str, "\t%.3g", nu_rev_score);
#endif

	asserta(!isnan(nu_rev_score));
	asserta(!isinf(nu_rev_score));
	asserta(!isnan(score));
	asserta(!isinf(score));

	const float selfw = m_params->m_self_w;
	if (selfw > 0)
		{
		float self_score = (m_self_rev_scores[i] + m_self_rev_scores[j])/2;
		score -= selfw*self_score;
		asserta(!isnan(m_self_rev_scores[i]));
		asserta(!isnan(m_self_rev_scores[j]));
		asserta(!isinf(score));
		asserta(!isnan(score));
		asserta(!isinf(score));
#if WRITE_TS_TERMS
		Psa(ts_str, "\t%.3g", self_score);
#endif
		}

	if (m_params->m_lddt_w > 0)
		{
		float lddt = flat_getlddt_muscle_some_floats(
			TD.m_pos_is, L_i,
			TD.m_pos_js, L_j,
			nmatch, distmx_i, distmx_j,
			TD.m_considered_vec, TD.m_preserved_vec);
		score += m_params->m_lddt_w*lddt*500;
		asserta(!isnan(lddt));
		asserta(!isinf(lddt));
		asserta(!isnan(score));
		asserta(!isinf(score));
#if WRITE_TS_TERMS
		Psa(ts_str, "\t%.3g", lddt);
#endif
		}

	if (m_params->m_lddtx_w > 0)
		{
		asserta(false);//lddtx screwed up?
		float L = (L_i + L_j)/2.0f + 50;
		float Lfactor = float(ncol)/L;

		float lddt = flat_getlddt_muscle_some_floats( // TODO this is not lddtx!?
			TD.m_pos_is, L_i,
			TD.m_pos_js, L_j,
			nmatch, distmx_i, distmx_j,
			TD.m_considered_vec, TD.m_preserved_vec);
		asserta(!isnan(lddt));
		asserta(!isinf(lddt));
		score += m_params->m_lddtx_w*lddt*500*Lfactor;
		asserta(!isnan(score));
		asserta(!isinf(score));
		}

	if (m_params->m_dali_w > 0)
		{
		float dali = flat_get_dali4(
			TD.m_pos_is, L_i,
			TD.m_pos_js, L_j,
			nmatch, distmx_i, distmx_j);
		asserta(!isnan(dali));
		asserta(!isinf(dali));
		score += m_params->m_dali_w*dali*10;
		asserta(!isnan(score));
		asserta(!isinf(score));
#if WRITE_TS_TERMS
		Psa(ts_str, "\t%.3g", dali);
#endif
		}

	if (m_params->m_dalix_w > 0)
		{
		Die("TODO");
		//float dalix = flat_get_dalix(
		//	label_i, label_j, path,
		//	lo_i, L_i, lo_j, L_j,
		//	distmx_i, distmx_j, TD.m_colscores);
		//score += m_params->m_dalix_w*dalix*10;
		}
#if WRITE_TS_TERMS
	Psa(ts_str, "\t%.3g", score);
	s_ts_lock.lock();
	fprintf(s_fts, "%s\n", ts_str.c_str());
	s_ts_lock.unlock();
#endif

	asserta(!isnan(score));
	asserta(!isinf(score));
	m_Scores[pairidx] = score;
	}

void flat_bench2::update_params(
	const vector<string> &names,
	const vector<float> &values)
	{
	if (m_single_feature)
		return;
	asserta(m_cdvec != 0);
	vector<string> alpha_names;
	vector<float> weights;
	vector<string> scalar_names;
	vector<float> scalar_values;
	flat_classify_params(
		names, values, alpha_names,
		weights, scalar_names, scalar_values);

	m_params->set_scalars(scalar_names, scalar_values);

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
	m_params->apply_weights(NameToWeight);
	chain_data::update_pssms(*m_params, m_cdvec, m_look->get_ndom());
	if (m_params->need_self())
		set_mega_self_rev_scores();
	if (m_params->need_nu_self())
		set_nu_self_rev_scores();
	}

//fprintf(f, "%s/%s", label_i.c_str(), sf_i.c_str());	0
//fprintf(f, "\t%s/%s", label_j.c_str(), sf_j.c_str());	1
//fprintf(f, "\t%u", L_i);	2
//fprintf(f, "\t%u", L_j);	3
//fprintf(f, "\t%u", lo_i);	4
//fprintf(f, "\t%u", lo_j);	5
//fprintf(f, "\t%u", lo_i_rev);	6
//fprintf(f, "\t%u", lo_j_rev);	7
//fprintf(f, "\t%s", cigar.c_str());	8
//fprintf(f, "\t%s", cigar_rev.c_str());	9
//fprintf(f, "\t%.1f", mega_fwd_score);	10
//fprintf(f, "\t%.1f", mega_score_rev);	11

void flat_bench2::load_mega_paths(const string &fn)
	{
	Alloc();
	const uint NQ = SIZE(m_Labels);
	const uint npair = triangle_get_K(NQ);
	for (uint i = 0; i < npair; ++i)
		m_Scores[i] = 0;

	m_mega_path_is.resize(npair);
	m_mega_path_js.resize(npair);
	m_mega_path_lo_i_fwds.resize(npair);
	m_mega_path_lo_j_fwds.resize(npair);
	m_mega_path_cigar_fwds.resize(npair);
	m_mega_path_score_fwds.resize(npair);
	m_mega_path_lo_i_revs.resize(npair);
	m_mega_path_lo_j_revs.resize(npair);
	m_mega_path_cigar_revs.resize(npair);
	m_mega_path_score_revs.resize(npair);

	FILE *f = OpenStdioFile(fn);
	string line;
	vector<string> flds;
	while (ReadLineStdioFile(f, line))
		{
		Split(line, flds, '\t');
		asserta(flds.size() == 12);
		const string &dom_i = flds[0];
		const string &dom_j = flds[1];
		uint i = m_look->get_domidx(dom_i);
		uint j = m_look->get_domidx(dom_j);
		uint k = triangle_ij_to_k2(i, j, NQ);

		m_mega_path_is[k] = i;
		m_mega_path_js[k] = j;
		m_mega_path_lo_i_fwds[k] = StrToUint(flds[4]);
		m_mega_path_lo_j_fwds[k] = StrToUint(flds[5]);
		m_mega_path_lo_i_revs[k] = StrToUint(flds[6]);
		m_mega_path_lo_j_revs[k] = StrToUint(flds[7]);
		m_mega_path_cigar_fwds[k] = flds[8];
		m_mega_path_cigar_revs[k] = flds[9];
		float score_fwd = StrToFloatf(flds[10]);
		m_mega_path_score_fwds[k] = score_fwd;
		m_mega_path_score_revs[k] = StrToFloatf(flds[11]);

		asserta(!isnan(score_fwd));
		asserta(!isinf(score_fwd));
		m_Scores[k] = float(score_fwd);
		}
	CloseStdioFile(f);
	}

void flat_bench2::load_single_feature(
	const string &fastafn,
	const string &logoddsfn,
	float gap_open, float gap_ext)
	{
	ProgressLog("load_single_feature(%s, %s)",
		fastafn.c_str(), logoddsfn.c_str());
	asserta(gap_open > 0);
	asserta(gap_ext >= 0);

	vector<float> logodds;
	uint alpha_size = flat_params::read_logodds(logoddsfn, logodds);
	m_feature_logodds = myalloc(float, alpha_size*alpha_size);
	for (uint i = 0; i < alpha_size*alpha_size; ++i)
		m_feature_logodds[i] = logodds[i];

	m_feature_gap_open = -gap_open;
	m_feature_gap_ext = -gap_ext;
	m_feature_alpha_size = alpha_size;
	m_single_feature = true;

	asserta(fastafn != "");
	asserta(logoddsfn != "");
	asserta(alpha_size > 1 && alpha_size <= 36);
	const uint8_t *char_to_letter =
		(alpha_size == 20 ? g_CharToLetterAmino : g_CharToLetterMu);

	SeqDB DB;
	DB.FromFasta(fastafn);
	DB.ToLetters(char_to_letter);
	DB.TruncLabels();
	DB.SetLabelToIndex();
	const uint ndom = m_look->get_ndom();
	m_feature_codeseq_vec = myalloc(uint8_t *, ndom);
	m_feature_codeseq_lengths = myalloc(uint, ndom);
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		const string &dom = m_look->get_dom(domidx);
		uint seqidx = DB.GetSeqIndex(dom);
		uint L = DB.GetSeqLength(seqidx);
		m_feature_codeseq_vec[domidx] = myalloc(uint8_t, L);
		memcpy(m_feature_codeseq_vec[domidx],
			(const uint8_t *) DB.m_Seqs[seqidx].c_str(), L);
		m_feature_codeseq_lengths[domidx] = L;
		}
	}

#if 0
void cmd_flat_bench2()
	{
#if WRITE_TS_TERMS
	s_fts = CreateStdioFile("ts.tmp");
#endif
	asserta(!optset_dope);
	asserta(!optset_subdope);

	const bool output_nu_paths = optset_output2;
	const bool input_mega_paths = optset_input2;
	if (output_nu_paths)
		s_f_nu_paths = CreateStdioFile(opt(output2));
	if (optset_output3)
		s_f_mega_paths = CreateStdioFile(opt(output3));

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
	flat_params params;
	params.set_scalars(scalar_names, scalar_values);
	params.init_from_alphadir(alphadir, alpha_names);

	vector<flat_chain_t *> chains;
	read_flat_chains(g_Arg1, chains);

	flat_bench2 FB;
	FB.m_params = &params;
	FB.ReadLookup(opt(lookup));
	if (optset_logodds)
		{
		asserta(optset_gapopen);
		asserta(optset_gapext);
		FB.load_single_feature(
			g_Arg1, opt(logodds),
			(float) opt(gapopen), (float) opt(gapext));
		}
	else
		{
		FB.load_chains(chains);
		FB.write_nu_hexfasta(opt(hexfasta));
		}

	FB.update_params(param_names, param_values);
	FB.m_nu_only = opt(nuonly);
	FB.m_timealn = opt(timealn);

	FB.m_params->logme();

	string varstr2;
	flat_make_varstr(params, varstr2);
	Log("varstr2=\n");
	Log("%s\n", varstr2.c_str());

	vector<string> peaker_spec_lines;
	flat_make_peaker_spec_const(params, peaker_spec_lines);
	Log("\n# const spec\n");
	for (uint i = 0; i < uint(peaker_spec_lines.size()); ++i)
		Log("%s\n", peaker_spec_lines[i].c_str());

	flat_make_peaker_spec_range(params, peaker_spec_lines);
	Log("\n# range spec\n");
	for (uint i = 0; i < uint(peaker_spec_lines.size()); ++i)
		Log("%s\n", peaker_spec_lines[i].c_str());

	uint nthread = GetRequestedThreadCount();
	thread_affinity ta;
	bool pin = opt(no_thread_pin) ? false : ta.shouldPin(nthread);
	if (input_mega_paths)
		{
		FB.m_input_mega_paths = true;
		FB.load_mega_paths(opt(input2));
		//FB.SetScoreOrder();
		//FB.Bench();
		FB.search(nthread, pin);
		FB.SetScoreOrder();
		FB.Bench();
		return;
		}

	FB.m_output_nu_paths = output_nu_paths;
	FB.search(nthread, pin);
	if (output_nu_paths)
		{
		CloseStdioFile(s_f_nu_paths);
		return;
		}
	CloseStdioFile(s_f_mega_paths);

	FB.SetScoreOrder();
	FB.Bench();
	FB.WriteHits(opt(output), opt(include_self), opt(triangle), opt(include_fam));

	double align_count = double(FB.m_aln_count);
	double mega_fwd_test_count = double(FB.m_mega_fwd_test_count);
	double mega_fwd_pass_count = double(FB.m_mega_fwd_pass_count);
	double nu_fwd_reject_count = double(FB.m_nu_fwd_reject_count);
	double nu_combined_reject_count = double(FB.m_nu_combined_reject_count);
	double mega_passed_pct = GetPct(mega_fwd_pass_count, mega_fwd_test_count);
	double ts_pct = GetPct(mega_fwd_pass_count, align_count);
	ProgressLog("Mu filter fwd %.1f%%, combined %.1f%%, total %.1f%%, passed %u\n",
		GetPct(nu_fwd_reject_count, align_count),
		GetPct(nu_combined_reject_count, align_count),
		GetPct(nu_fwd_reject_count+nu_combined_reject_count, align_count),
		FB.m_nu_pass_count.load());
	ProgressLog("Mega fwd filter passed %.1f%%\n", mega_passed_pct);
#if WRITE_TS_TERMS
	CloseStdioFile(s_fts);
#endif
	}
#endif
