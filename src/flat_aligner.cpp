#include "myutils.h"
#include "flat_params.h"
#include "flat_features.h"
#include "flat_aligner.h"
#include "flat_helpers.h"
#include "flat_alignx.h"
#include "cigar.h"

void flat_aligner::alloc()
	{
	const uint nfeat = flat_features::m_nfeat;

	m_pssmT = myalloc(float, m_maxL*flat_features::get_sum_alpha_sizes());
	m_pssm_reverseT = myalloc(float, m_maxL*flat_features::get_sum_alpha_sizes());

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
	fill_flat_pssm_reversed(profT, LT, flat_features::m_nfeat,
		flat_features::m_alpha_sizes,
		flat_features::m_feature_block_offsets,
		flat_features::m_weighted_logoddsvec,
		m_pssmT);
	}

void flat_aligner::cache_reverseT(const string &labelT, const uint8_t *profT, uint LT)
	{
	asserta(LT < m_maxL);
	m_labelT = labelT;
	m_profT = profT;
	m_LT = LT;
	fill_flat_pssm_reversed(profT, LT, flat_features::m_nfeat,
		flat_features::m_alpha_sizes,
		flat_features::m_feature_block_offsets,
		flat_features::m_weighted_logoddsvec,
		m_pssm_reverseT);
	}

void flat_aligner::cacheT(const string &labelT, const uint8_t *profT, uint LT)
	{
	asserta(LT < m_maxL);
	m_labelT = labelT;
	m_profT = profT;
	m_LT = LT;
	fill_flat_pssm(profT, LT, flat_features::m_nfeat,
		flat_features::m_alpha_sizes,
		flat_features::m_feature_block_offsets,
		flat_features::m_weighted_logoddsvec,
		m_pssmT);
	}

void flat_aligner::alignQ(const string &labelQ, const uint8_t *profQ, uint LQ)
	{
	m_labelQ = labelQ;
	m_profQ = profQ;
	m_LQ = LQ;
	m_score = sw_flat_pssm(
		m_scratch_rows, m_TB, m_scratch_pssms,
		profQ, LQ,
		m_pssmT, m_LT, 
		flat_features::m_feature_block_offsets, flat_features::m_nfeat,
		-flat_params::m_open, 
		-flat_params::m_ext,
		m_loQ, m_loT, m_path_buffer, m_ncol);
	m_reverse_score_set = false;
	}

void flat_aligner::align_reverse()
	{
	assert(m_profQ != 0);
	assert(m_LQ > 0);
	m_reverse_score = sw_flat_pssm(
		m_scratch_rows, m_TB, m_scratch_pssms,
		m_profQ, m_LQ,
		m_pssm_reverseT, m_LT, 
		flat_features::m_feature_block_offsets,
		flat_features::m_nfeat,
		-flat_params::m_open, 
		-flat_params::m_ext,
		m_loQ, m_loT, m_path_buffer, m_ncol);
	m_reverse_score_set = true;
	}

float flat_aligner::get_self_rev_score(
	const string &labelQ, const uint8_t *profQ, uint LQ)
	{
	cacheT_reversed(labelQ + ".rev", profQ, LQ);
	alignQ(labelQ, profQ, LQ);
	m_reverse_score_set = false;
	return m_score;
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
	uint nfeat = flat_features::m_nfeat;
	assert(nfeat > 0);
	const uint32_t *alpha_sizes = flat_features::m_alpha_sizes;
	const vector<string> &feature_names = flat_features::m_feature_names;
	const vector<string> &symbolsvec = flat_features::m_symbolsvec;

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
		fprintf(f, "  %8.8s*%2u", feature_names[fi].c_str(), alpha_sizes[fi]);
		fprintf(f, "  %s\n", m_labelQ.c_str());

		fprintf(f, "%s\n", annot_row.c_str());

		fprintf(f, "%s", feature_rowT.c_str());
		fprintf(f, "  %8.8s*%2u", feature_names[fi].c_str(), alpha_sizes[fi]);
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
