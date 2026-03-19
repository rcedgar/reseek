#pragma once
#include "alpha.h"

float sw_flat_pssm(
	float *__restrict scratch_rows,
	uint8_t *__restrict TB,
	const float ** __restrict scratch_ppsms,
	const uint8_t *__restrict profA, uint LA,
	const float *__restrict pssm, uint LB,
	const uint32_t * __restrict feature_block_offsets,
	uint nfeat,
	float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
	string &Path);

void profiles2faprof(
	const string &fn,
	const vector<string> &feature_names,
	const vector<uint> &alpha_sizes,
	const vector<string> &labels,
	const vector<vector<uint8_t> > &profiles);

void read_profiles_faprof(
	const string &fn,
	vector<string> &feature_names,
	vector<uint> &alpha_sizes,
	vector<string> &labels,
	vector<vector<uint8_t> > &profiles);

void read_profiles_and_logoddsmxvec(
	const string &specfn,
	vector<string> &feature_names,
	vector<uint> &alpha_sizes,
	vector<string> &labels,
	vector<vector<uint8_t> > &profiles,
	vector<vector<float> > &logoddsmxvec);

void check_profiles(
	vector<vector<uint8_t> > &profiles,
	vector<uint> &alpha_sizes);

uint32_t get_flat_pssm_feature_block_offsets(
	const uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	uint32_t * __restrict feature_block_offsets);

void fill_flat_pssm(
	const uint8_t * __restrict profQ,
	uint32_t LQ,
	uint32_t nfeat,
	const uint32_t * __restrict alpha_sizes,
	const uint32_t * __restrict feature_block_offsets,
	const float *const * __restrict weighted_logoddsmxvec,
	float * __restrict pssm);

void fill_smx_using_flat_pssm(
	const uint8_t * __restrict profA,
	uint32_t LA,
	uint32_t LB,
	uint32_t nfeat,
	const uint32_t * __restrict feature_block_offsets,
	const float * __restrict pssm,
	float * __restrict smx);

static inline const uint8_t *get_letter2char(uint alpha_size)
	{
	return (alpha_size == 20 ? g_LetterToCharAmino : g_LetterToCharMu);
	}

static inline const uint8_t *get_char2letter(uint alpha_size)
	{
	return (alpha_size == 20 ? g_CharToLetterAmino : g_CharToLetterMu);
	}
