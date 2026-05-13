#include "myutils.h"
#include "cigar.h"
#include "flat_nu_aligner.h"
#include "final_nu_matrix.h"

parasail_matrix_t flat_nu_aligner::m_matrix;
int flat_nu_aligner::m_open = INT_MAX;	// penalty > 0
int flat_nu_aligner::m_ext = INT_MAX;	// penalty > 0
int flat_nu_aligner::m_saturated_score = INT_MAX;

atomic<uint32_t> flat_nu_aligner::m_aln_count;
atomic<uint32_t> flat_nu_aligner::m_saturated_count;

/***
$src/reseek_tune2/bash/final_nu_matrix.bash
                                 vvvvvvvvvvvvvv--- scale pre-built into matrix
intopen=2.90E+01;intext=3.00E+00;scale=8.81E+00;aa4=5.15E-01;pm2=2.84E-01;sec32=2.00E-01;
***/
bool flat_nu_aligner::init()
	{
	m_open = 29;
	m_ext = 3;
	int minscore = 0;
	int maxscore = 0;
	for (uint i = 0; i < 256*256; ++i)
		{
		int score = s_final_nu_matrix[i];
		if (i == 0 || score < minscore) minscore = score;
		if (i == 0 || score > maxscore) maxscore = score;
		}
	m_matrix.size = 256;
	m_matrix.length = 256;
	m_matrix.type = PARASAIL_MATRIX_TYPE_SQUARE;
	m_matrix.matrix = s_final_nu_matrix;
	m_matrix.min = minscore;
	m_matrix.max = maxscore;
#if DEBUG
	int *mapper = myalloc(int, 256);
	memset(mapper, 0, 256*sizeof(int));
	for (int i = 0; i < 256; ++i)
		mapper[i] = i;
	m_matrix.mapper = mapper;
#else
	m_matrix.mapper = 0;
#endif
	return true;
	}
static bool init_done = flat_nu_aligner::init();

void flat_nu_aligner::set_query(
	const string &labelQ,
	const byte *codeseqQ,
	uint LQ)
	{
	m_labelQ = labelQ;
	m_LQ = LQ;
	m_codeseqQ = codeseqQ;
	if (m_parasail_profQ != 0)
		parasail_profile_free(m_parasail_profQ);
	m_parasail_profQ = parasail_profile_create_avx_256_16(
		(const char *) codeseqQ, LQ, &m_matrix);
	}

int flat_nu_aligner::align_score_only(
	const string &labelT,
	const byte *codeseqT,
	uint LT)
	{
	++m_aln_count;
	clear_aln();
	m_labelT = labelT;
	m_codeseqT = codeseqT;
	m_LT = LT;
	if (m_result != 0)
		parasail_result_free(m_result);

	m_result = parasail_sw_striped_profile_avx2_256_16(
		m_parasail_profQ, (const char *) codeseqT, LT, m_open, m_ext);
	if (m_result->flag & PARASAIL_FLAG_SATURATED)
		{
		m_score = m_saturated_score;
		++m_saturated_count;
		}
	else
		m_score = m_result->score;
	return m_score;
	}

int flat_nu_aligner::align_path(
	const string &labelT,
	const byte *codeseqT,
	uint LT)
	{
	clear_aln();
	m_labelT = labelT;
	m_codeseqT = codeseqT;
	m_LT = LT;
	if (m_result != 0)
		parasail_result_free(m_result);

	m_result = parasail_sw_trace_striped_profile_avx2_256_8(
		m_parasail_profQ,(const char *) codeseqT, LT, m_open, m_ext);
	asserta(!(m_result->flag & PARASAIL_FLAG_SATURATED));
	m_score = m_result->score;

	parasail_cigar_t* cig = parasail_result_get_cigar_extra(
		m_result,
		(const char *) m_codeseqQ, m_LQ,
		(const char *) m_codeseqT, m_LT,
		&m_matrix, 1, 0);

	char *cig_str = parasail_cigar_decode(cig);
	m_loQ = (uint) cig->beg_query;
	m_loT = (uint) cig->beg_ref;
	ExpandParaCigar_reverseDI(cig_str, m_path);
	free(cig_str);
	parasail_cigar_free(cig);

#if DEBUG
	{
	assert(m_loQ == 0);
	assert(m_loT == 0);
	uint M, D, I;
	GetPathCounts(m_path, M, D, I);
	asserta(m_path.back() == 'M');
	asserta(M + D <= m_LQ);
	asserta(M + I <= m_LT);
	}
#endif
	return m_score;
	}
