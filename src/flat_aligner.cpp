#include "myutils.h"
#include "flat_params.h"
#include "flat_alphas.h"
#include "flat_aligner.h"
#include "flat_helpers.h"
#include "flat_alignx.h"
#include "paralign.h"
#include "cigar.h"

//atomic<uint> flat_aligner::m_nu_filter_reject_count;
atomic<uint> flat_aligner::m_aln_count;

void flat_aligner::alloc()
	{
	const uint nfeat = flat_alphas::m_nfeat;

	m_pssmT = myalloc(float, m_maxL*flat_alphas::get_sum_alpha_sizes());
	m_pssm_reverseT = myalloc(float, m_maxL*flat_alphas::get_sum_alpha_sizes());

	m_scratch_rows = myalloc(float, 2*m_maxL + 2);
	m_scratch_pssms = myalloc(const float *, nfeat);
	m_TB = myalloc(uint8_t, m_maxL*m_maxL);
	m_path_buffer = myalloc(char, 2*m_maxL);
	}

void flat_aligner::freemem()
	{
	myfree(m_pssmT);		m_pssmT = 0;
	myfree(m_pssm_reverseT);m_pssm_reverseT = 0;
	myfree(m_scratch_rows);	m_scratch_rows = 0;
	myfree(m_scratch_pssms);m_scratch_pssms = 0;
	myfree(m_TB);			m_TB = 0;
	myfree(m_path_buffer);	m_path_buffer = 0;
	}

void flat_aligner::cacheT_reversed(const string &labelT, const uint8_t *profT, uint LT)
	{
	asserta(LT < m_maxL);
	m_labelT = labelT;
	m_profT = profT;
	m_LT = LT;
	fill_flat_pssm_reversed(profT, LT, flat_alphas::m_nfeat,
		flat_alphas::m_alpha_sizes,
		flat_alphas::m_feature_block_offsets,
		flat_alphas::m_weighted_logoddsvec,
		m_pssmT);
	}

void flat_aligner::cache_reverseT(
	const string &labelT,
	const uint8_t *profT,
	const uint8_t *nu_codeseq_rev,
	uint LT)
	{
	asserta(LT < m_maxL);
	m_labelT = labelT;
	m_profT = profT;
	m_LT = LT;
	fill_flat_pssm_reversed(profT, LT, flat_alphas::m_nfeat,
		flat_alphas::m_alpha_sizes,
		flat_alphas::m_feature_block_offsets,
		flat_alphas::m_weighted_logoddsvec,
		m_pssm_reverseT);
	//if (m_nu_filter)
	//	{
	//	assert(m_pa);
	//	assert(nu_codeseq_rev);
	//	m_pa->SetQueryProfile_rev(nu_codeseq_rev, LT);
	//	}
	}

void flat_aligner::cacheT(
	const string &labelT,
	const uint8_t *profT,
	const uint8_t *nu_codeseq,
	uint LT)
	{
	asserta(LT < m_maxL);
	m_labelT = labelT;
	m_profT = profT;
	m_LT = LT;

	//if (m_nu_only || m_nu_filter)
	//	{
	//	//TODO query<->target
	//	assert(m_pa);
	//	assert(nu_codeseq);
	//	m_pa->SetQueryProfile(m_labelT, nu_codeseq, LT);
	//	if (m_nu_only)
	//		return;
	//	}

	fill_flat_pssm(profT, LT, flat_alphas::m_nfeat,
		flat_alphas::m_alpha_sizes,
		flat_alphas::m_feature_block_offsets,
		flat_alphas::m_weighted_logoddsvec,
		m_pssmT);

	//if (m_nu_filter)
	//	{
	//	assert(m_pa);
	//	assert(nu_codeseq);
	//	m_pa->SetQueryProfile(m_labelT, nu_codeseq, LT);
	//	}
	}

void flat_aligner::alignQ(
	const string &labelQ,
	const uint8_t *profQ,
	const uint8_t *nu_codeseqQ,
	uint LQ)
	{
	clear_align();

	m_labelQ = labelQ;
	m_profQ = profQ;
	m_LQ = LQ;

	++m_aln_count;
	//if (m_nu_only || m_nu_filter)
	//	{
	//	assert(m_pa);
	//	assert(nu_codeseqQ);
	//	m_pa->Align_ScoreOnly(labelQ, nu_codeseqQ, LQ);
	//	m_score = float(m_pa->m_Score);
	//	if (m_nu_only)
	//		return;
	//	if (m_pa->m_Score < flat_params::m_min_nu_fwd_score)
	//		{
	//		m_nu_filter_reject = true;
	//		++m_nu_filter_reject_count;
	//		return;
	//		}
	//	}
	//m_nu_filter_reject = false;
	m_score = sw_flat_pssm(
		m_scratch_rows, m_TB, m_scratch_pssms,
		profQ, LQ,
		m_pssmT, m_LT, 
		flat_alphas::m_feature_block_offsets,
		flat_alphas::m_nfeat,
		-flat_params::m_open, 
		-flat_params::m_ext,
		m_loQ, m_loT, m_path_buffer, m_ncol);
#if DEBUG
	validate_path();
#endif
	}

void flat_aligner::align_reverse()
	{
	assert(m_profQ != 0);
	assert(m_LQ > 0);
	m_reverse_score = sw_flat_pssm(
		m_scratch_rows, m_TB, m_scratch_pssms,
		m_profQ, m_LQ,
		m_pssm_reverseT, m_LT, 
		flat_alphas::m_feature_block_offsets,
		flat_alphas::m_nfeat,
		-flat_params::m_open, 
		-flat_params::m_ext,
		m_loQ, m_loT, m_path_buffer, m_ncol);
	m_reverse_score_set = true;
	}

float flat_aligner::get_self_rev_score(
	const string &labelQ, const uint8_t *profQ, uint LQ)
	{
	m_reverse_score_set = false;
	const uint nfeat = flat_alphas::get_nfeat();
	uint8_t *revprofQ = myalloc(uint8_t, LQ*nfeat);
	flat_reverse_profile(profQ, LQ, nfeat, revprofQ);
	cacheT(labelQ + ".rev", revprofQ, 0, LQ);
	alignQ(labelQ, profQ, 0, LQ);
	myfree(revprofQ);
	return m_score;
	}

void flat_aligner::validate_path() const
	{
	asserta(m_loQ < m_LQ);
	asserta(m_loT < m_LT);
	uint posQ = m_loQ;
	uint posT = m_loT;
	for (uint i = 0; i < m_ncol; ++i)
		{
		char c = m_path_buffer[i];
		if (c == 'M')
			{
			asserta(posQ < m_LQ);
			asserta(posT < m_LT);
			++posQ;
			++posT;
			}
		else if (c == 'D')
			{
			asserta(posQ < m_LQ);
			++posQ;
			}
		else if (c == 'I')
			{
			asserta(posT < m_LT);
			++posT;
			}
		else
			asserta(false);
		}
	}

uint flat_aligner::get_match_count() const
	{
	uint m = 0;
	for (uint i = 0; i < m_ncol; ++i)
		if (m_path_buffer[i] == 'M') ++m;
	return m;
	}

uint flat_aligner::get_path_str(string &path) const
	{
	uint m = 0;
	path.clear();
	path.reserve(m_ncol);
	for (uint i = 0; i < m_ncol; ++i)
		{
		char c = m_path_buffer[i];
		path.push_back(c);
		if (c == 'M') ++m;
		}
	return m;
	}

void flat_aligner::write_tsv(FILE *f) const
	{
	if (f == 0) return;
	string CIGAR;
	PathToCIGAR(m_path_buffer, CIGAR);

	fprintf(f, "%s", m_labelQ.c_str());
	fprintf(f, "\t%s", m_labelT.c_str());
	fprintf(f, "\t%.4g", m_score);
	fprintf(f, "\t%u", m_loQ);
	fprintf(f, "\t%u", m_LQ);
	fprintf(f, "\t%u", m_loT);
	fprintf(f, "\t%u", m_LT);
	fprintf(f, "\t%s", CIGAR.c_str());
	fprintf(f, "\n");
	}

void flat_aligner::write_aln(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, "\n");
	uint nfeat = flat_alphas::m_nfeat;
	assert(nfeat > 0);
	const uint32_t *alpha_sizes = flat_alphas::m_alpha_sizes;
	const vector<string> &alpha_names = flat_alphas::m_alpha_names;
	const vector<string> &symbolsvec = flat_alphas::m_symbolsvec;

	vector<string> feature_rowsQ(nfeat);
	vector<string> feature_rowsT(nfeat);
	uint hiQ = 0;
	uint hiT = 0;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		string feature_rowQ = feature_rowsQ[fi];
		string feature_rowT = feature_rowsT[fi];
		string annot_row;

		feature_rowQ.reserve(m_ncol);
		feature_rowT.reserve(m_ncol);
		annot_row.reserve(m_ncol);
		uint alpha_size = alpha_sizes[fi];
		uint posQ = m_loQ;
		uint posT = m_loT;
		const uint8_t *letter2char = get_letter2char(alpha_size);
		for (uint col = 0; col < m_ncol; ++col)
			{
			char c = m_path_buffer[col];
			if (c == 'M')
				{
				assert(posQ < m_LQ);
				assert(posT < m_LT);
				uint8_t codeQ = m_profQ[fi*m_LQ + posQ];
				uint8_t codeT = m_profT[fi*m_LT + posT];
				annot_row += (codeQ == codeT) ? '|' :
					symbolsvec[fi][codeQ*alpha_size + codeT];
				}
			else
				annot_row += ' ';

			if (c == 'M' || c == 'D')
				{
				assert(posQ < m_LQ);
				uint8_t codeQ = m_profQ[fi*m_LQ + posQ];
				feature_rowQ += letter2char[codeQ];
				++posQ;
				}
			else
				feature_rowQ += '-';

			if (c == 'M' || c == 'I')
				{
				assert(posT < m_LT);
				uint8_t codeT = m_profT[fi*m_LT + posT];
				feature_rowT += letter2char[codeT];
				++posT;
				}
			else
				feature_rowT += '-';
			}
		fprintf(f, "\n");
		fprintf(f, "%s", feature_rowQ.c_str());
		fprintf(f, "  %8.8s*%2u", alpha_names[fi].c_str(), alpha_sizes[fi]);
		fprintf(f, "  %s\n", m_labelQ.c_str());

		fprintf(f, "%s\n", annot_row.c_str());

		fprintf(f, "%s", feature_rowT.c_str());
		fprintf(f, "  %8.8s*%2u", alpha_names[fi].c_str(), alpha_sizes[fi]);
		fprintf(f, "  %s\n", m_labelT.c_str());

		if (fi == 0)
			{
			hiQ = posQ;
			hiT = posT;
			}
		}
	fprintf(f, "score %.1f", m_score);
	fprintf(f, ", Q %u-%u(%u)", m_loQ+1, hiQ, m_LQ);
	fprintf(f, ", T %u-%u(%u)", m_loT+1, hiT, m_LT);
	fprintf(f, "\n");
	}
