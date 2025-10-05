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

static const int parasail_mu_map[256] = {
0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};

// reseek -train_feature ../traindata/train.fa2 -db ../data/scop40.bca -feature Mu -alpha_size 36 -wildcard no
// reseek.exe -musubstmx Mu.36 -output musubstmx.txt -log musubstmx.log [0cd02b5]
static const int parasail_mu_[36*36] = {
 3,-3, 0, 1,-6,-2, 2,-5,-2, 2,-3, 0, 0,-5,-2, 1,-4,-2, 0,-3,-2,-1,-5,-4, 0,-4,-3,-1,-2,-2,-2,-4,-3,-1,-3,-3,  // 0
-3, 4, 2,-4, 1, 0,-3, 1, 0,-4, 3, 1,-4, 0, 0,-4, 1, 0,-3, 2, 0,-5, 0,-2,-4, 0,-1,-3, 1,-1,-4,-2,-2,-2,-1,-2,  // 1
 0, 2, 3,-1,-2, 1, 0, 0, 1,-1, 1, 2,-2,-2, 0,-1,-1, 1,-1, 0, 1,-3,-3,-1,-2,-1, 0,-2,-1,-1,-3,-3,-2,-2,-2,-2,  // 2
 1,-4,-1, 3,-5,-1, 2,-5,-1, 0,-6,-2, 2,-5,-1, 1,-5,-2,-1,-4,-3, 1,-6,-2, 0,-5,-3,-2,-4,-3,-1,-4,-3,-1,-4,-3,  // 3
-6, 1,-2,-5, 2, 0,-6, 1,-1,-6, 0,-3,-6, 1,-1,-6, 0,-2,-5,-1,-4,-6, 0,-2,-5,-1,-3,-4,-2,-3,-5,-1,-2,-4,-2,-3,  // 4
-2, 0, 1,-1, 0, 2,-2, 0, 1,-3,-1, 0,-2, 0, 1,-2, 0, 0,-4,-2,-2,-3,-1, 0,-3,-1,-1,-4,-2,-2,-4,-2,-1,-3,-3,-2,  // 5
 2,-3, 0, 2,-6,-2, 2,-4,-1, 1,-4,-1, 1,-6,-2, 1,-4,-1, 0,-4,-2, 0,-5,-3, 0,-4,-2,-1,-3,-3,-1,-4,-3,-1,-4,-3,  // 6
-5, 1, 0,-5, 1, 0,-4, 2, 1,-5, 1,-2,-5, 0,-1,-5, 1,-1,-5, 0,-3,-5,-1,-2,-5, 0,-2,-4,-1,-3,-4,-2,-2,-3,-1,-2,  // 7
-2, 0, 1,-1,-1, 1,-1, 1, 2,-3,-1, 0,-2,-1, 0,-2, 0, 1,-3,-1,-1,-3,-2,-1,-2,-1, 0,-3,-2,-2,-4,-3,-2,-3,-2,-1,  // 8
 2,-4,-1, 0,-6,-3, 1,-5,-3, 2,-3, 0, 1,-5,-2, 1,-4,-2, 1,-3,-1, 0,-5,-3, 0,-4,-2, 0,-3,-2,-1,-4,-3,-1,-3,-3,  // 9
-3, 3, 1,-6, 0,-1,-4, 1,-1,-3, 4, 1,-4, 1, 0,-4, 2, 0,-3, 3, 1,-4, 0,-1,-3, 1, 0,-2, 1,-1,-3,-1,-2,-2, 0,-2,  // 10
 0, 1, 2,-2,-3, 0,-1,-2, 0, 0, 1, 2,-2,-2, 0,-1, 0, 1, 0, 1, 2,-2,-2, 0,-1,-1, 1,-2,-1, 0,-3,-3,-2,-2,-2,-1,  // 11
 0,-4,-2, 2,-6,-2, 1,-5,-2, 1,-4,-2, 3,-5, 0, 1,-5,-1, 0,-4,-2, 2,-5,-1, 1,-4,-2,-1,-3,-3, 0,-4,-2,-1,-4,-3,  // 12
-5, 0,-2,-5, 1, 0,-6, 0,-1,-5, 1,-2,-5, 3, 1,-5, 2,-1,-4, 0,-3,-5, 2, 0,-5, 1,-2,-3,-1,-3,-5, 0,-2,-4,-1,-3,  // 13
-2, 0, 0,-1,-1, 1,-2,-1, 0,-2, 0, 0, 0, 1, 2,-2, 0, 1,-3,-1, 0,-1, 0, 1,-2, 0, 0,-3,-2,-1,-2,-2, 0,-3,-2,-1,  // 14
 1,-4,-1, 1,-6,-2, 1,-5,-2, 1,-4,-1, 1,-5,-2, 2,-4,-1, 1,-3,-2, 1,-5,-2, 1,-4,-1,-1,-3,-2,-1,-4,-3, 0,-3,-2,  // 15
-4, 1,-1,-5, 0, 0,-4, 1, 0,-4, 2, 0,-5, 2, 0,-4, 3, 1,-4, 1,-1,-4, 1, 0,-3, 2, 0,-3, 0,-2,-4,-1,-2,-3, 0,-1,  // 16
-2, 0, 1,-2,-2, 0,-1,-1, 1,-2, 0, 1,-1,-1, 1,-1, 1, 2,-2, 0, 0,-2,-1, 0,-1, 0, 1,-3,-1,-1,-3,-3,-1,-2,-2,-1,  // 17
 0,-3,-1,-1,-5,-4, 0,-5,-3, 1,-3, 0, 0,-4,-3, 1,-4,-2, 2,-2, 0, 1,-4,-2, 1,-4,-1, 1,-2,-1,-1,-3,-2, 0,-3,-2,  // 18
-3, 2, 0,-4,-1,-2,-4, 0,-1,-3, 3, 1,-4, 0,-1,-3, 1, 0,-2, 4, 1,-4, 1, 0,-3, 2, 0,-2, 3, 0,-4,-1,-1,-2, 1,-1,  // 19
-2, 0, 1,-3,-4,-2,-2,-3,-1,-1, 1, 2,-2,-3, 0,-2,-1, 0, 0, 1, 2,-1,-2, 0,-1, 0, 1,-1, 0, 1,-2,-3,-1,-1,-1, 0,  // 20
-1,-5,-3, 1,-6,-3, 0,-5,-3, 0,-4,-2, 2,-5,-1, 1,-4,-2, 1,-4,-1, 3,-4, 0, 2,-4,-1, 0,-2,-2, 1,-3,-1, 0,-3,-2,  // 21
-5, 0,-3,-6, 0,-1,-5,-1,-2,-5, 0,-2,-5, 2, 0,-5, 1,-1,-4, 1,-2,-4, 3, 1,-4, 2,-1,-4, 0,-2,-4, 1,-1,-3, 0,-2,  // 22
-4,-2,-1,-2,-2, 0,-3,-2,-1,-3,-1, 0,-1, 0, 1,-2, 0, 0,-2, 0, 0, 0, 1, 2,-1, 1, 1,-2,-1,-1,-1,-1, 1,-2,-1, 0,  // 23
 0,-4,-2, 0,-5,-3, 0,-5,-2, 0,-3,-1, 1,-5,-2, 1,-3,-1, 1,-3,-1, 2,-4,-1, 2,-3, 0, 0,-2,-1, 0,-3,-2, 1,-3,-1,  // 24
-4, 0,-1,-5,-1,-1,-4, 0,-1,-4, 1,-1,-4, 1, 0,-4, 2, 0,-4, 2, 0,-4, 2, 1,-3, 3, 1,-3, 1,-1,-3, 0,-1,-3, 1, 0,  // 25
-3,-1, 0,-3,-3,-1,-2,-2, 0,-2, 0, 1,-2,-2, 0,-1, 0, 1,-1, 0, 1,-1,-1, 1, 0, 1, 2,-2,-1, 0,-2,-2, 0,-1,-1, 0,  // 26
-1,-3,-2,-2,-4,-4,-1,-4,-3, 0,-2,-2,-1,-3,-3,-1,-3,-3, 1,-2,-1, 0,-4,-2, 0,-3,-2, 2, 0, 1, 1,-2,-1, 1,-2,-1,  // 27
-2, 1,-1,-4,-2,-2,-3,-1,-2,-3, 1,-1,-3,-1,-2,-3, 0,-1,-2, 3, 0,-2, 0,-1,-2, 1,-1, 0, 4, 2,-1, 1, 0,-1, 2, 1,  // 28
-2,-1,-1,-3,-3,-2,-3,-3,-2,-2,-1, 0,-3,-3,-1,-2,-2,-1,-1, 0, 1,-2,-2,-1,-1,-1, 0, 1, 2, 2,-1,-1, 0, 0, 0, 1,  // 29
-2,-4,-3,-1,-5,-4,-1,-4,-4,-1,-3,-3, 0,-5,-2,-1,-4,-3,-1,-4,-2, 1,-4,-1, 0,-3,-2, 1,-1,-1, 3,-2, 0, 2,-2, 0,  // 30
-4,-2,-3,-4,-1,-2,-4,-2,-3,-4,-1,-3,-4, 0,-2,-4,-1,-3,-3,-1,-3,-3, 1,-1,-3, 0,-2,-2, 1,-1,-2, 3, 1,-2, 2, 0,  // 31
-3,-2,-2,-3,-2,-1,-3,-2,-2,-3,-2,-2,-2,-2, 0,-3,-2,-1,-2,-1,-1,-1,-1, 1,-2,-1, 0,-1, 0, 0, 0, 1, 2, 0, 1, 1,  // 32
-1,-2,-2,-1,-4,-3,-1,-3,-3,-1,-2,-2,-1,-4,-3, 0,-3,-2, 0,-2,-1, 0,-3,-2, 1,-3,-1, 1,-1, 0, 2,-2, 0, 2,-1, 0,  // 33
-3,-1,-2,-4,-2,-3,-4,-1,-2,-3, 0,-2,-4,-1,-2,-3, 0,-2,-3, 1,-1,-3, 0,-1,-3, 1,-1,-2, 2, 0,-2, 2, 1,-1, 3, 1,  // 34
-3,-2,-2,-3,-3,-2,-3,-2,-1,-3,-2,-1,-3,-3,-1,-2,-1,-1,-2,-1, 0,-2,-2, 0,-1, 0, 0,-1, 1, 1, 0, 0, 1, 0, 1, 2,  // 35
};

parasail_matrix_t parasail_mu_matrix = {
	"mu",  // name
	parasail_mu_,  // matrix
	parasail_mu_map,  // mapper
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
		  result, SeqA, LA, SeqB, LB, &parasail_mu_matrix, 1, 0);

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
		++m_ParasailSaturateCount;
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
	const uint LA = SIZE(*m_MuLettersA);
	m_ProfPara = parasail_profile_create_avx_256_8(A, LA, &parasail_mu_matrix);

	m_MuRevA.clear();
	m_MuRevA.reserve(LA);
	for (uint i = 0; i < LA; ++i)
		m_MuRevA.push_back((*m_MuLettersA)[LA-i-1]);
	const char *AR = (const char *) m_MuRevA.data();
	m_ProfParaRev = parasail_profile_create_avx_256_8(AR, LA, &parasail_mu_matrix);
	EndTimer(SetMuQP_Para);
	}

float DSSAligner::AlignMuParaBags(const ChainBag &BagA, const ChainBag &BagB)
	{
	asserta(BagA.m_ptrProfPara != 0);
	asserta(BagA.m_ptrProfParaRev != 0);
	asserta(BagB.m_ptrMuLetters != 0);

	uint LB = BagB.m_ptrChain->GetSeqLength();
	asserta(SIZE(*BagB.m_ptrMuLetters) == LB);

	const int Open = m_Params->m_ParaMuGapOpen;
	const int Ext = m_Params->m_ParaMuGapExt;
	const float OmegaFwd = m_Params->m_OmegaFwd;

	const char *B = (const char *) BagB.m_ptrMuLetters->data();

	const parasail_profile_t * const restrict profile =
	  (const parasail_profile_t * const restrict) BagA.m_ptrProfPara;
	parasail_result_t* result =
	  parasail_sw_striped_profile_avx2_256_8(profile, B, LB, Open, Ext);
	if (result->flag & PARASAIL_FLAG_SATURATED)
		{
		++m_ParasailSaturateCount;
		result->score = 777;
		}
	float fwd_score = (float) result->score;
	if (fwd_score < OmegaFwd)
		{
		parasail_result_free(result);
		return 0;
		}

	const parasail_profile_t * const restrict profile_rev =
	  (const parasail_profile_t * const restrict) BagA.m_ptrProfParaRev;
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
