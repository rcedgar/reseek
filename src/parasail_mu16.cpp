#include "myutils.h"
#include "parasail.h"
#include "dssaligner.h"
#include "cigar.h"
#include "timing.h"

static const int parasail_mu_map[256] = {
0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,};

/***
C:\src\2025-10_reseek_tune\2025-11-28_retrain_para_mu_sweep\report.tsv
1.263	mints10.maxts60.scale5	open=20	ext=11
1.274	evaluerange.scale2		open=11	ext=4
1.290	mints10.maxts25.scale2	open=8	ext=5 <<<

[5f0310b]
reseek \
	-hjmumx hjmumx.spec \
	-db ../data/scop40.bca \
	-input2 ../out/a.bca \
	-lookup ../out/a.lookup \
	-fixmubyteseq \
	-scale 2 \
	-traintps ../big_fa2/tp.mints10.maxts25.fa2 -log mints10.maxts25.scale2

                                     vvvvvvvvvvvvvv
Paralign::LogMatrix() open 0, ext 0, min -19, max 8, size 36, length 36
***/
int parasail_mu_16[36*36] = {
  5, -6,  1,  3,-13, -5,  4,-10, -3,  2,-10, -2,  0,-17, -6,  1,-12, -5, -1,-10, -5, -3,  0, -8, -2,-15, -6, -3,-13, -6, -5,  0,-10, -4,  0, -9,  // 0
 -6,  8,  4, -9,  2,  0, -6,  3,  1,-10,  5,  0,-13,  0, -3,-10,  1, -2,-10,  3, -2,  0,  0, -4,  0,  0, -4,-14,  0, -5,  0, -4, -7,  0, -3, -6,  // 1
  1,  4,  6, -1, -3,  2,  0,  0,  3, -2,  0,  2, -4, -5, -1, -3, -2,  0, -4, -1,  1, -6, -6, -2, -4, -3, -1, -6, -4, -2, -9, -8, -5, -7, -6, -4,  // 2
  3, -9, -1,  7, -9, -1,  5, -9, -2,  0,-13, -5,  4,-11, -2,  2,-11, -4, -3,  0, -7,  2,-11, -4, -1,-12, -5, -4,  0, -8, -1,  0, -7, -3,  0, -9,  // 3
-13,  2, -3, -9,  6,  2,-11,  4,  0,-17, -1, -7,-13,  2, -2,-15,  0, -5,  0, -3, -9,-15,  0, -4,-17, -2, -7,-19, -5,-12,-17, -1, -6,  0, -3, -9,  // 4
 -5,  0,  2, -1,  2,  6, -2,  2,  4, -7, -3, -2, -3, -1,  2, -5, -1,  0, -9, -5, -4, -5, -2,  0, -6, -3, -2,-12, -7, -6, -8, -4, -2, -9, -5, -4,  // 5
  4, -6,  0,  5,-11, -2,  5, -7,  0,  1, -9, -3,  2,-13, -4,  3,-10, -3, -1,-13, -5,  0,-15, -5,  0,-11, -4, -4,-13, -8, -3,  0, -7, -2,-14, -7,  // 6
-10,  3,  0, -9,  4,  2, -7,  6,  2,-13,  0, -4,-13,  0, -2,-11,  2, -2,-16, -1, -6,-12,  0, -4,-12,  0, -3,  0, -4, -8,  0, -3, -6,-16, -2, -6,  // 7
 -3,  1,  3, -2,  0,  4,  0,  2,  5, -6, -3,  0, -5, -3,  0, -3, -1,  1, -7, -4, -2, -6, -4, -1, -5, -2, -1,-10, -6, -5, -8, -6, -4, -7, -4, -3,  // 8
  2,-10, -2,  0,-17, -7,  1,-13, -6,  3, -6, -1,  0,-14, -6,  2,-10, -4,  2, -9, -3, -1,-18, -7,  0,-12, -6,  0,-11, -5, -3,  0, -9, -2,-15, -8,  // 9
-10,  5,  0,-13, -1, -3, -9,  0, -3, -6,  6,  1,-10,  1, -2, -8,  2, -1, -8,  5,  0,-12,  0, -3, -9,  1, -2,-10,  3, -2,-12, -1, -4,-10,  0, -4,  // 10
 -2,  0,  2, -5, -7, -2, -3, -4,  0, -1,  1,  3, -4, -5, -1, -2, -2,  1, -2,  0,  2, -5, -6, -2, -3, -3,  0, -4, -2,  0, -6, -7, -3, -5, -4, -2,  // 11
  0,-13, -4,  4,-13, -3,  2,-13, -5,  0,-10, -4,  6, -8, -1,  3, -9, -3, -1,  0, -6,  4,-10, -2,  1,-11, -4, -3,  0, -8,  3,-13, -4,  0,-12, -6,  // 12
-17,  0, -5,-11,  2, -1,-13,  0, -3,-14,  1, -5, -8,  5,  1,-10,  3, -2,-18, -1, -6,-10,  3, -1,-14,  1, -4,-18, -2, -8,-10,  3, -2,-13,  0, -5,  // 13
 -6, -3, -1, -2, -2,  2, -4, -2,  0, -6, -2, -1, -1,  1,  4, -3,  0,  1, -7, -3, -2, -2,  0,  2, -4, -2,  0, -8, -5, -4, -3, -1,  1, -5, -2, -1,  // 14
  1,-10, -3,  2,-15, -5,  3,-11, -3,  2, -8, -2,  3,-10, -3,  3, -7, -2,  0,-11, -4,  1,-13, -4,  2,-10, -3, -2,-16, -6,  0,-15, -6,  0,-11, -5,  // 15
-12,  1, -2,-11,  0, -1,-10,  2, -1,-10,  2, -2, -9,  3,  0, -7,  4,  0,-12,  1, -4,-11,  1, -2,-10,  2, -2,-14,  0, -5,-10,  1, -2,-10,  2, -2,  // 16
 -5, -2,  0, -4, -5,  0, -3, -2,  1, -4, -1,  1, -3, -2,  1, -2,  0,  2, -5, -2, -1, -4, -3,  0, -3, -1,  1, -7, -4, -2, -5, -4, -1, -4, -2,  0,  // 17
 -1,-10, -4, -3,  0, -9, -1,-16, -7,  2, -8, -2, -1,-18, -7,  0,-12, -5,  4, -5,  0,  2,-12, -4,  3,-10, -3,  2, -8, -2,  0,  0, -6,  1,-11, -5,  // 18
-10,  3, -1,  0, -3, -5,-13, -1, -4, -9,  5,  0,  0, -1, -3,-11,  1, -2, -5,  8,  3, -9,  3,  0, -7,  4,  0, -7,  5,  0, -9, -1, -3, -8,  1, -2,  // 19
 -5, -2,  1, -7, -9, -4, -5, -6, -2, -3,  0,  2, -6, -6, -2, -4, -4, -1,  0,  3,  4, -3, -3,  1, -1,  0,  2, -2,  1,  2, -5, -5, -2, -3, -3,  0,  // 20
 -3,  0, -6,  2,-15, -5,  0,-12, -6, -1,-12, -5,  4,-10, -2,  1,-11, -4,  2, -9, -3,  6, -5,  1,  4, -8, -1,  0,-12, -5,  4,-10, -1,  2,-10, -3,  // 21
  0,  0, -6,-11,  0, -2,-15,  0, -4,-18,  0, -6,-10,  3,  0,-13,  1, -3,-12,  3, -3, -5,  7,  3,-10,  4, -1,-14,  0, -5, -8,  4,  0,-11,  1, -3,  // 22
 -8, -4, -2, -4, -4,  0, -5, -4, -1, -7, -3, -2, -2, -1,  2, -4, -2,  0, -4,  0,  1,  1,  3,  5, -1,  1,  2, -6, -3, -2, -1,  0,  3, -3, -1,  0,  // 23
 -2,  0, -4, -1,-17, -6,  0,-12, -5,  0, -9, -3,  1,-14, -4,  2,-10, -3,  3, -7, -1,  4,-10, -1,  5, -7,  0,  1,-11, -4,  2,-13, -3,  3, -8, -2,  // 24
-15,  0, -3,-12, -2, -3,-11,  0, -2,-12,  1, -3,-11,  1, -2,-10,  2, -1,-10,  4,  0, -8,  4,  1, -7,  6,  1,-12,  2, -3, -9,  2, -1, -8,  3, -1,  // 25
 -6, -4, -1, -5, -7, -2, -4, -3, -1, -6, -2,  0, -4, -4,  0, -3, -2,  1, -3,  0,  2, -1, -1,  2,  0,  1,  3, -5, -2,  0, -3, -3,  0, -2, -1,  1,  // 26
 -3,-14, -6, -4,-19,-12, -4,  0,-10,  0,-10, -4, -3,-18, -8, -2,-14, -7,  2, -7, -2,  0,-14, -6,  1,-12, -5,  4, -5,  0,  2,-13, -4,  2, -9, -3,  // 27
-13,  0, -4,  0, -5, -7,-13, -4, -6,-11,  3, -2,  0, -2, -5,-16,  0, -4, -8,  5,  1,-12,  0, -3,-11,  2, -2, -5,  7,  2, -9,  2, -1, -7,  3,  0,  // 28
 -6, -5, -2, -8,-12, -6, -8, -8, -5, -5, -2,  0, -8, -8, -4, -6, -5, -2, -2,  0,  2, -5, -5, -2, -4, -3,  0,  0,  2,  4, -2, -3,  1, -1, -1,  2,  // 29
 -5,  0, -9, -1,-17, -8, -3,  0, -8, -3,-12, -6,  3,-10, -3,  0,-10, -5,  0, -9, -5,  4, -8, -1,  2, -9, -3,  2, -9, -2,  7, -6,  1,  4, -7, -1,  // 30
  0, -4, -8,  0, -1, -4,  0, -3, -6,  0, -1, -7,-13,  3, -1,-15,  1, -4,  0, -1, -5,-10,  4,  0,-13,  2, -3,-13,  2, -3, -6,  7,  2, -9,  4,  0,  // 31
-10, -7, -5, -7, -6, -2, -7, -6, -4, -9, -4, -3, -4, -2,  1, -6, -2, -1, -6, -3, -2, -1,  0,  3, -3, -1,  0, -4, -1,  1,  1,  2,  5, -1,  1,  3,  // 32
 -4,  0, -7, -3,  0, -9, -2,-16, -7, -2,-10, -5,  0,-13, -5,  0,-10, -4,  1, -8, -3,  2,-11, -3,  3, -8, -2,  2, -7, -1,  4, -9, -1,  4, -6,  0,  // 33
  0, -3, -6,  0, -3, -5,-14, -2, -4,-15,  0, -4,-12,  0, -2,-11,  2, -2,-11,  1, -3,-10,  1, -1, -8,  3, -1, -9,  3, -1, -7,  4,  1, -6,  6,  2,  // 34
 -9, -6, -4, -9, -9, -4, -7, -6, -3, -8, -4, -2, -6, -5, -1, -5, -2,  0, -5, -2,  0, -3, -3,  0, -2, -1,  1, -3,  0,  2, -1,  0,  3,  0,  2,  4,  // 35
};

parasail_matrix_t parasail_mu_matrix = {
	"mu",  // name
	parasail_mu_16,  // matrix
	parasail_mu_map,  // mapper
	36,  // size
	8, // max
	-19,  // min
	NULL, // user_matrix
	PARASAIL_MATRIX_TYPE_SQUARE,  // type
	36,  // length
	"ABCDEFGHIJKLMNOPQRSTUVWZYZabcdefghij", // alphabet
	NULL // query
};

void Log_parasail_mu_matrix(const parasail_matrix_t &mx)
	{
	Log("name   %s\n", mx.name);
	Log("size   %d\n", mx.size);
	Log("max    %d\n", mx.max);
	Log("min    %d\n", mx.min);
	Log("length %d\n", mx.length);
	Log("type   %d\n", mx.type);
	Log("user_mx  %p\n", mx.user_matrix);
	Log("mapper\n");
	for (uint i = 0; i < 36; ++i)
		Log(" %d", mx.mapper[i]);
	Log("\n");
	Log("matrix\n");
	for (uint i = 0; i < 36; ++i)
		{
		for (uint j = 0; j < 36; ++j)
			Log(" %3d", mx.matrix[36*i + j]);
		Log("\n");
		}
	}

void Log_parasail_profile(const parasail_profile_t &prof)
	{
	Log("Log_parasail_profile()\n");
	Log("s1Len %d\n", prof.s1Len);
	Log("s1 ");
	for (int i = 0; i < prof.s1Len; ++i)
		Log(" %d", prof.s1[i]);
	Log("\n");
	const int8_t *ptr_score = (const int8_t *) prof.profile8.score;
	int k = 0;
	for (int i = 0; i < 36; ++i)
		{
		Log("%3d  |%3d| ", i, prof.s1[i]);
		for (uint j = 0; j < (uint) prof.s1Len; ++j)
			Log(" %3d", ptr_score[k++]);
		Log("\n");
		}
	}

float DSSAligner::AlignMuQP_Para()
	{
	m_MuFwdScore = 0;
	m_MuRevScore = 0;
	StartTimer(SWPara);
	uint LA = SIZE(*m_MuLettersA);
	uint LB = SIZE(*m_MuLettersB);
	const int Open = DSSParams::m_ParaMuGapOpen;
	const int Ext = DSSParams::m_ParaMuGapExt;
	const float OmegaFwd = DSSParams::m_OmegaFwd;

	const char *B = (const char *) m_MuLettersB->data();

	const parasail_profile_t * const restrict profile =
	  (const parasail_profile_t * const restrict) m_ProfPara;
	parasail_result_t* result =
	  parasail_sw_striped_profile_avx2_256_16(profile, B, LB, Open, Ext);
#if 0
	{
	Log_parasail_mu_matrix(parasail_mu_matrix);
	const byte *m_Q = m_MuLettersA->data();
	const byte *m_T = m_MuLettersB->data();
	uint m_LQ = LA;
	uint m_LT = LB;
	Log("QL %u, TL %u\n", m_LQ, m_LT);
	Log("Q: ");
	for (uint i = 0; i < m_LQ; ++i)
		Log(" %u", m_Q[i]);
	Log("\n");
	Log("T: ");
	for (uint i = 0; i < m_LT; ++i)
		Log(" %u", m_T[i]);
	Log("\n");
	Log("score %d\n", result->score);
	}
#endif

	if (result->flag & PARASAIL_FLAG_SATURATED)
		{
		++m_ParasailSaturateCount;
		result->score = 777;
		}
	m_MuFwdScore = (float) result->score;
	if (m_MuFwdScore < OmegaFwd)
		{
		parasail_result_free(result);
		EndTimer(SWPara);
		return 0;
		}

	const parasail_profile_t * const restrict profile_rev =
	  (const parasail_profile_t * const restrict) m_ProfParaRev;
	parasail_result_t* result_rev =
	  parasail_sw_striped_profile_avx2_256_16(profile_rev, B, LB, Open, Ext);
	m_MuRevScore = (float) result_rev->score;

	EndTimer(SWPara);
	if (result_rev->flag & PARASAIL_FLAG_SATURATED)
		result_rev->score = 777;
	m_MuFwdMinusRevScore = m_MuFwdScore - m_MuRevScore;
	parasail_result_free(result);
	parasail_result_free(result_rev);
	return m_MuFwdMinusRevScore;
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
	m_ProfPara = parasail_profile_create_avx_256_16(A, LA, &parasail_mu_matrix);
	const parasail_profile_t *prof =
		(parasail_profile_t *) m_ProfPara;
#if 0
	{
	Log_parasail_profile(*prof);
	Log("ProfQ:");
	const byte *m_Q = m_MuLettersA->data();
	uint m_LQ = LA;
	Log("QL %u\n", m_LQ);
	Log("Q: ");
	for (uint i = 0; i < m_LQ; ++i)
		Log(" %u", m_Q[i]);
	Log("\n");
	const byte *ptrProf = (const byte *) m_ProfPara;
	for (uint i = 0; i < sizeof(parasail_profile_t); ++i)
		Log(" %02x", ptrProf[i]);
	Log("\n");
	}
#endif
	m_MuRevA.clear();
	m_MuRevA.reserve(LA);
	for (uint i = 0; i < LA; ++i)
		m_MuRevA.push_back((*m_MuLettersA)[LA-i-1]);
	const char *AR = (const char *) m_MuRevA.data();
	m_ProfParaRev = parasail_profile_create_avx_256_16(AR, LA, &parasail_mu_matrix);
	EndTimer(SetMuQP_Para);
	}

float DSSAligner::AlignMuParaBags(const ChainBag &BagA, const ChainBag &BagB)
	{
	asserta(BagA.m_ptrProfPara != 0);
	asserta(BagA.m_ptrProfParaRev != 0);
	asserta(BagB.m_ptrMuLetters != 0);

	uint LB = BagB.m_ptrChain->GetSeqLength();
	asserta(SIZE(*BagB.m_ptrMuLetters) == LB);

	const int Open = DSSParams::m_ParaMuGapOpen;
	const int Ext = DSSParams::m_ParaMuGapExt;
	const float OmegaFwd = DSSParams::m_OmegaFwd;

	const char *B = (const char *) BagB.m_ptrMuLetters->data();

	const parasail_profile_t * const restrict profile =
	  (const parasail_profile_t * const restrict) BagA.m_ptrProfPara;
	parasail_result_t* result =
	  parasail_sw_striped_profile_avx2_256_16(profile, B, LB, Open, Ext);
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
	  parasail_sw_striped_profile_avx2_256_16(profile_rev, B, LB, Open, Ext);
	float rev_score = (float) result_rev->score;

	EndTimer(SWPara);
	if (result_rev->flag & PARASAIL_FLAG_SATURATED)
		result_rev->score = 777;
	float Score = fwd_score - rev_score;
	parasail_result_free(result);
	parasail_result_free(result_rev);
	return Score;
	}
