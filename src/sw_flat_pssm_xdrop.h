#pragma once

#include "myutils.h"
#include "tracebit.h"

// Forward X-drop Smith–Waterman on flat multi-feature PSSM (0-based coords).
// Explores i >= posQ, j >= posT. Open/Ext must be <= 0 (same as sw_flat_pssm).
// TB length LQ*LT; scratch_rows length 2*LT+3; scratch_ppsms length nfeat.
float sw_flat_pssm_xdrop_fwd(
	float *__restrict scratch_rows,
	uint8_t *__restrict TB,
	const float ** __restrict scratch_ppsms,
	const uint8_t *__restrict profQ, uint LQ,
	const float *__restrict pssmT, uint LT,
	const uint32_t *__restrict feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext,
	uint &loQ, uint &loT,
	char *path_buffer, uint &ncol);

float sw_flat_pssm_xdrop_fwd_scoreonly(
	float *__restrict scratch_rows,
	const float ** __restrict scratch_ppsms,
	const uint8_t *__restrict profQ, uint LQ,
	const float *__restrict pssmT, uint LT,
	const uint32_t *__restrict feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext);

// Reverse extension from (posQ-1, posT-1); score 0 if posQ==0 or posT==0.
float sw_flat_pssm_xdrop_bwd(
	float *__restrict scratch_rows,
	uint8_t *__restrict TB,
	const float ** __restrict scratch_ppsms,
	const uint8_t *__restrict profQ, uint LQ,
	const float *__restrict pssmT, uint LT,
	const uint32_t *__restrict feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext,
	uint &loQ, uint &loT,
	char *path_buffer, uint &ncol);

float sw_flat_pssm_xdrop_bwd_scoreonly(
	float *__restrict scratch_rows,
	const float ** __restrict scratch_ppsms,
	const uint8_t *__restrict profQ, uint LQ,
	const float *__restrict pssmT, uint LT,
	const uint32_t *__restrict feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext);

// Fwd from (posQ,posT) + bwd from (posQ-1,posT-1), merged path (XDropHSP-style).
float sw_flat_pssm_xdrop_hsp(
	float *__restrict scratch_rows,
	uint8_t *__restrict TB_fwd,
	uint8_t *__restrict TB_bwd,
	const float ** __restrict scratch_ppsms,
	const uint8_t *__restrict profQ, uint LQ,
	const float *__restrict pssmT, uint LT,
	const uint32_t *__restrict feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext,
	uint &loQ, uint &loT,
	char *path_buffer, uint &ncol);

// Reference: same X-drop rule, no j-band (validates banded implementation).
float sw_flat_pssm_xdrop_fwd_ref(
	float *__restrict scratch_rows,
	uint8_t *__restrict TB,
	const float ** __restrict scratch_ppsms,
	const uint8_t *__restrict profQ, uint LQ,
	const float *__restrict pssmT, uint LT,
	const uint32_t *__restrict feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float Open, float Ext,
	uint8_t *__restrict active, // optional, LQ*LT, 1 = visited
	uint &loQ, uint &loT,
	char *path_buffer, uint &ncol);

// --- tests (call from cmd_test_sw_flat_pssm_xdrop or STANDALONE driver) ---
void test_sw_flat_pssm_xdrop_banded_vs_ref(
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	float open, float ext,
	uint ntrial);

void test_sw_flat_pssm_xdrop_X_huge_vs_full(
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	float open, float ext);

void test_sw_flat_pssm_xdrop_active_paths(
	const uint8_t *profQ, uint LQ,
	const float *pssmT, uint LT,
	const uint32_t *feature_block_offsets,
	uint nfeat,
	uint posQ, uint posT,
	float X, float open, float ext);
