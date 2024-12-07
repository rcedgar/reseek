#include "myutils.h"
#include "parasail.h"
#include "dssaligner.h"
#include "cigar.h"
#include "timing.h"

/***
This source code is partly derived from the parasail library
https://github.com/jeffdaily/parasail
Author jeffrey.daily@gmail.com
Copyright (c) 2015 Battelle Memorial Institute.
***/

static const int parasail_combo_map[256] = {
0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};

static const int parasail_combo_[36*36] = {
 3,-1, 1, 1,-3,-1, 2,-2, 0, 1,-3,-1,-1,-5,-3, 0,-4,-2, 0,-4,-2,-2,-6,-4,-1,-5,-3,-1,-5,-3,-3,-7,-5,-2,-6,-4,  // 0
-1, 4, 2,-3, 1,-1,-2, 2, 0,-3, 2, 0,-5, 0,-2,-4, 0,-2,-4, 1,-1,-6,-1,-3,-5,-1,-3,-5, 0,-2,-7,-2,-4,-6,-2,-4,  // 1
 1, 2, 3,-1,-1, 0, 0, 0, 1,-1, 0, 1,-3,-2,-1,-2,-2,-1,-2,-1, 0,-4,-3,-2,-3,-3,-2,-3,-2,-1,-5,-4,-3,-4,-4,-3,  // 2
 1,-3,-1, 3, 0, 1, 2,-2, 0,-1,-5,-3, 2,-2,-1, 0,-4,-2,-2,-6,-4, 0,-3,-2,-1,-5,-3,-3,-7,-5,-1,-4,-3,-2,-6,-4,  // 3
-3, 1,-1, 0, 4, 2,-2, 3, 1,-5, 0,-2,-2, 2, 0,-4, 1,-1,-6,-1,-3,-3, 1,-1,-5, 0,-2,-7,-2,-4,-4, 0,-2,-6,-1,-3,  // 4
-1,-1, 0, 1, 2, 3, 0, 1, 2,-3,-2,-1,-1, 0, 1,-2,-1, 0,-4,-3,-2,-2,-1, 0,-3,-2,-1,-5,-4,-3,-3,-2,-1,-4,-3,-2,  // 5
 2,-2, 0, 2,-2, 0, 3,-1, 1, 0,-4,-2, 0,-4,-2, 1,-3,-1,-1,-5,-3,-1,-5,-3, 0,-4,-2,-2,-6,-4,-2,-6,-4,-1,-5,-3,  // 6
-2, 2, 0,-2, 3, 1,-1, 3, 1,-4, 0,-2,-4, 1,-1,-3, 1,-1,-5,-1,-3,-5, 0,-2,-4, 0,-2,-6,-2,-4,-6,-1,-3,-5,-1,-3,  // 7
 0, 0, 1, 0, 1, 2, 1, 1, 2,-2,-2,-1,-2,-1, 0,-1,-1, 0,-3,-3,-2,-3,-2,-1,-2,-2,-1,-4,-4,-3,-4,-3,-2,-3,-3,-2,  // 8
 1,-3,-1,-1,-5,-3, 0,-4,-2, 2,-2, 0, 0,-4,-2, 1,-3,-1, 1,-3,-1,-1,-5,-3, 0,-4,-2, 1,-3,-1,-2,-5,-4,-1,-5,-3,  // 9
-3, 2, 0,-5, 0,-2,-4, 0,-2,-2, 3, 1,-4, 0,-2,-3, 1,-1,-3, 2, 0,-5, 0,-2,-4, 0,-2,-3, 1,-1,-5,-1,-3,-5, 0,-2,  // 10
-1, 0, 1,-3,-2,-1,-2,-2,-1, 0, 1, 2,-2,-2,-1,-1,-1, 0,-1, 0, 1,-3,-2,-1,-2,-2,-1,-1,-1, 0,-4,-3,-2,-3,-2,-1,  // 11
-1,-5,-3, 2,-2,-1, 0,-4,-2, 0,-4,-2, 2,-1, 0, 1,-3,-1,-1,-5,-3, 1,-2,-1, 0,-4,-2,-2,-5,-4, 1,-3,-1, 0,-4,-2,  // 12
-5, 0,-2,-2, 2, 0,-4, 1,-1,-4, 0,-2,-1, 3, 1,-3, 2, 0,-5, 0,-2,-2, 2, 0,-4, 1,-1,-5,-1,-3,-3, 1,-1,-4, 0,-2,  // 13
-3,-2,-1,-1, 0, 1,-2,-1, 0,-2,-2,-1, 0, 1, 2,-1, 0, 1,-3,-2,-1,-1, 0, 1,-2,-1, 0,-4,-3,-2,-1,-1, 0,-2,-2,-1,  // 14
 0,-4,-2, 0,-4,-2, 1,-3,-1, 1,-3,-1, 1,-3,-1, 2,-2, 0, 0,-4,-2, 0,-4,-2, 1,-3,-1,-1,-5,-3, 0,-4,-2, 0,-3,-2,  // 15
-4, 0,-2,-4, 1,-1,-3, 1,-1,-3, 1,-1,-3, 2, 0,-2, 2, 0,-4, 0,-2,-4, 1,-1,-3, 1,-1,-5, 0,-2,-4, 0,-2,-3, 1,-1,  // 16
-2,-2,-1,-2,-1, 0,-1,-1, 0,-1,-1, 0,-1, 0, 1, 0, 0, 1,-2,-2,-1,-2,-1, 0,-1,-1, 0,-3,-2,-1,-2,-2,-1,-2,-1, 0,  // 17
 0,-4,-2,-2,-6,-4,-1,-5,-3, 1,-3,-1,-1,-5,-3, 0,-4,-2, 3,-1, 1, 0,-3,-2, 1,-2,-1, 1,-2,-1,-1,-5,-3, 0,-4,-2,  // 18
-4, 1,-1,-6,-1,-3,-5,-1,-3,-3, 2, 0,-5, 0,-2,-4, 0,-2,-1, 3, 1,-3, 1,-1,-2, 2, 0,-2, 2, 0,-5, 0,-2,-4, 0,-2,  // 19
-2,-1, 0,-4,-3,-2,-3,-3,-2,-1, 0, 1,-3,-2,-1,-2,-2,-1, 1, 1, 2,-2,-1, 0,-1, 0, 1,-1, 0, 1,-3,-2,-1,-2,-2,-1,  // 20
-2,-6,-4, 0,-3,-2,-1,-5,-3,-1,-5,-3, 1,-2,-1, 0,-4,-2, 0,-3,-2, 3,-1, 1, 2,-2, 0,-1,-5,-3, 2,-2, 0, 0,-4,-2,  // 21
-6,-1,-3,-3, 1,-1,-5, 0,-2,-5, 0,-2,-2, 2, 0,-4, 1,-1,-3, 1,-1,-1, 4, 2,-2, 2, 0,-5, 0,-2,-2, 2, 0,-4, 1,-1,  // 22
-4,-3,-2,-2,-1, 0,-3,-2,-1,-3,-2,-1,-1, 0, 1,-2,-1, 0,-2,-1, 0, 1, 2, 3, 0, 0, 1,-3,-2,-1, 0, 0, 1,-2,-1, 0,  // 23
-1,-5,-3,-1,-5,-3, 0,-4,-2, 0,-4,-2, 0,-4,-2, 1,-3,-1, 1,-2,-1, 2,-2, 0, 2,-1, 0, 0,-4,-2, 0,-4,-2, 1,-3,-1,  // 24
-5,-1,-3,-5, 0,-2,-4, 0,-2,-4, 0,-2,-4, 1,-1,-3, 1,-1,-2, 2, 0,-2, 2, 0,-1, 3, 1,-4, 0,-2,-4, 1,-1,-3, 2, 0,  // 25
-3,-3,-2,-3,-2,-1,-2,-2,-1,-2,-2,-1,-2,-1, 0,-1,-1, 0,-1, 0, 1, 0, 0, 1, 0, 1, 2,-2,-2,-1,-2,-1, 0,-1, 0, 1,  // 26
-1,-5,-3,-3,-7,-5,-2,-6,-4, 1,-3,-1,-2,-5,-4,-1,-5,-3, 1,-2,-1,-1,-5,-3, 0,-4,-2, 3,-1, 1, 0,-3,-2, 1,-3,-1,  // 27
-5, 0,-2,-7,-2,-4,-6,-2,-4,-3, 1,-1,-5,-1,-3,-5, 0,-2,-2, 2, 0,-5, 0,-2,-4, 0,-2,-1, 3, 1,-3, 1,-1,-3, 2, 0,  // 28
-3,-2,-1,-5,-4,-3,-4,-4,-3,-1,-1, 0,-4,-3,-2,-3,-2,-1,-1, 0, 1,-3,-2,-1,-2,-2,-1, 1, 1, 2,-2,-1, 0,-1, 0, 1,  // 29
-3,-7,-5,-1,-4,-3,-2,-6,-4,-2,-5,-4, 1,-3,-1, 0,-4,-2,-1,-5,-3, 2,-2, 0, 0,-4,-2, 0,-3,-2, 3,-1, 1, 2,-2,-1,  // 30
-7,-2,-4,-4, 0,-2,-6,-1,-3,-5,-1,-3,-3, 1,-1,-4, 0,-2,-5, 0,-2,-2, 2, 0,-4, 1,-1,-3, 1,-1,-1, 3, 1,-2, 2, 0,  // 31
-5,-4,-3,-3,-2,-1,-4,-3,-2,-4,-3,-2,-1,-1, 0,-2,-2,-1,-3,-2,-1, 0, 0, 1,-2,-1, 0,-2,-1, 0, 1, 1, 2,-1, 0, 1,  // 32
-2,-6,-4,-2,-6,-4,-1,-5,-3,-1,-5,-3, 0,-4,-2, 0,-3,-2, 0,-4,-2, 0,-4,-2, 1,-3,-1, 1,-3,-1, 2,-2,-1, 2,-2, 0,  // 33
-6,-2,-4,-6,-1,-3,-5,-1,-3,-5, 0,-2,-4, 0,-2,-3, 1,-1,-4, 0,-2,-4, 1,-1,-3, 2, 0,-3, 2, 0,-2, 2, 0,-2, 3, 1,  // 34
-4,-4,-3,-4,-3,-2,-3,-3,-2,-3,-2,-1,-2,-2,-1,-2,-1, 0,-2,-2,-1,-2,-1, 0,-1, 0, 1,-1, 0, 1,-1, 0, 1, 0, 1, 2,  // 35
};

parasail_matrix_t parasail_combo_matrix = {
	"combo",  // name
	parasail_combo_,  // matrix
	parasail_combo_map,  // mapper
	36,  // size
	4, // max
	-7,  // min
	NULL, // user_matrix
	PARASAIL_MATRIX_TYPE_SQUARE,  // type
	36,  // length
	"ABCDEFGHIJKLMNOPQRSTUVWZYZabcdefghij", // alphabet
	NULL // query
};

float DSSAligner::AlignMuQP_Para_Path(uint &LoA, uint &LoB, string &Path)
	{
	Path.clear();
	LoA = UINT_MAX;
	LoB = UINT_MAX;

	uint LA = SIZE(*m_MuLettersA);
	uint LB = SIZE(*m_MuLettersB);
	const int Open = m_Params->m_ParaMuGapOpen;
	const int Ext = m_Params->m_ParaMuGapExt;

	const char *B = (const char *) m_MuLettersB->data();

	const parasail_profile_t * const restrict profile =
	  (const parasail_profile_t * const restrict) m_ProfPara;

	StartTimer(SWPara);
	parasail_result_t* result =
	  parasail_sw_trace_striped_profile_avx2_256_8(profile, B, LB, Open, Ext);
	EndTimer(SWPara);

	float Score = (float) result->score;
	if (result->flag & PARASAIL_FLAG_SATURATED)
		Score = -1;
	else
		{
		const char *SeqA = (const char *) m_MuLettersA->data();
		const char *SeqB = (const char *) m_MuLettersB->data();
		parasail_cigar_t* cig = parasail_result_get_cigar_extra(
		  result, SeqA, LA, SeqB, LB, &parasail_combo_matrix, 1, 0);

		char *cig_str = parasail_cigar_decode(cig);
		ExpandParaCigar_reverseDI(cig_str, Path);

		LoA = (uint) cig->beg_query;
		LoB = (uint) cig->beg_ref;

		free(cig_str);
		parasail_cigar_free(cig);
		}
	parasail_result_free(result);
	return Score;
	}

float DSSAligner::AlignMuQP_Para()
	{
	StartTimer(SWPara);
	uint LA = SIZE(*m_MuLettersA);
	uint LB = SIZE(*m_MuLettersB);
	const int Open = m_Params->m_ParaMuGapOpen;
	const int Ext = m_Params->m_ParaMuGapExt;
	const float OmegaFwd = m_Params->m_OmegaFwd;

	const char *B = (const char *) m_MuLettersB->data();

	const parasail_profile_t * const restrict profile =
	  (const parasail_profile_t * const restrict) m_ProfPara;
	parasail_result_t* result =
	  parasail_sw_striped_profile_avx2_256_8(profile, B, LB, Open, Ext);
	if (result->flag & PARASAIL_FLAG_SATURATED)
		{
		m_StatsLock.lock();
		++m_ParasailSaturateCount;
		m_StatsLock.unlock();
		result->score = 777;
		}
	float fwd_score = (float) result->score;
	if (fwd_score < OmegaFwd)
		{
		parasail_result_free(result);
		EndTimer(SWPara);
		return 0;
		}

	const parasail_profile_t * const restrict profile_rev =
	  (const parasail_profile_t * const restrict) m_ProfParaRev;
	parasail_result_t* result_rev =
	  parasail_sw_striped_profile_avx2_256_8(profile_rev, B, LB, Open, Ext);
	float rev_score = (float) result_rev->score;

	EndTimer(SWPara);
	if (result_rev->flag & PARASAIL_FLAG_SATURATED)
		result_rev->score = 777;
	float Score = fwd_score - rev_score;
	parasail_result_free(result);
	parasail_result_free(result_rev);
	return Score;
	}

void DSSAligner::SetMuQP_Para()
	{
	StartTimer(SetMuQP_Para);
	if (m_ProfPara != 0)
		parasail_profile_free((parasail_profile_t *) m_ProfPara);
	if (m_ProfParaRev != 0)
		parasail_profile_free((parasail_profile_t *) m_ProfParaRev);
	const char *A = (const char *) m_MuLettersA->data();
	uint LA = SIZE(*m_MuLettersA);
	m_ProfPara = parasail_profile_create_avx_256_8(A, LA, &parasail_combo_matrix);

	vector<byte> ARev = *m_MuLettersA;
	const char *AR = (const char *) ARev.data();
	reverse(ARev.begin(), ARev.end());
	m_ProfParaRev = parasail_profile_create_avx_256_8(AR, LA, &parasail_combo_matrix);
	EndTimer(SetMuQP_Para);
	}
