#include "myutils.h"
#include "dss.h"
#include "chaq.h"
#include "flat_base.h"
#include "flat_chain.h"
#include "pdbchain.h"
#include "flat_distmx.h"
#include "pdbfilescanner.h"
#include "flat_chain_reader.h"

static const uint MAXL = 1000;

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

float SWFast(XDPMem &Mem, const float * const *SMxData, uint LA, uint LB,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path);

static uint s_nfeat;
static uint *s_alpha_sizes;
static uint *s_feature_block_offsets;
static float *s_pssm_i;
static float *s_smx;
static float **s_weighted_logoddsmxvec;
static uint s_L_i;
static string s_label_i;

static void cache_i_flat(const string &label, const uint8_t *prof_i, uint L_i)
	{
	s_label_i = label;
	s_L_i = L_i;
	fill_flat_pssm(prof_i, L_i, s_nfeat, s_alpha_sizes,
		s_feature_block_offsets, s_weighted_logoddsmxvec, s_pssm_i);
	}

static float align_j_flat(const string &label, const uint8_t *prof_j, uint L_j)
	{
	fill_smx_using_flat_pssm(prof_j, L_j, s_L_i, s_nfeat,
		s_feature_block_offsets, s_pssm_i, s_smx);
	return 0;
	}

static void cache_i_old(const string &label, const uint8_t *prof_i, uint L_i)
	{
	}

static float align_j_old(const string &label, const uint8_t *prof_j, uint L_j)
	{
	return 0;
	}

void cmd_test_sw()
	{
	const string &specfn = g_Arg1;
	vector<string> labels;
	vector<vector<uint8_t> > profiles;
	vector<string> feature_names;
	vector<uint> alpha_sizes;
	vector<vector<float> > logoddsmxvec;
	read_profiles_and_logoddsmxvec(
		specfn,
		feature_names,
		alpha_sizes,
		labels,
		profiles,
		logoddsmxvec);

	s_nfeat = SIZE(feature_names);
	const uint nprof = SIZE(labels);

	asserta(SIZE(alpha_sizes) == s_nfeat);
	asserta(SIZE(profiles) == nprof);

	s_smx = myalloc(float, MAXL*MAXL);
	s_alpha_sizes = alpha_sizes.data();

	check_profiles(profiles, alpha_sizes);

	s_weighted_logoddsmxvec = myalloc(float *, s_nfeat);
	for (uint fi = 0; fi < s_nfeat; ++fi)
		{
		uint AS = alpha_sizes[fi];
		asserta(AS >= 2 && AS < 256);
		s_weighted_logoddsmxvec[fi] = logoddsmxvec[fi].data();
		}

	s_feature_block_offsets = myalloc(uint32_t, s_nfeat);
	const uint32_t rows_per_pos =
		get_flat_pssm_feature_block_offsets(s_nfeat,
			s_alpha_sizes, s_feature_block_offsets);

	for (uint i = 0; i < nprof; ++i)
		{
		uint L_i = SIZE(profiles[i]);
		if (L_i > MAXL) continue;
		cache_i_flat(labels[i], profiles[i].data(), L_i);
		cache_i_old(labels[i], profiles[i].data(), L_i);

		for (uint j = 0; j < nprof; ++j)
			{
			uint L_j = SIZE(profiles[j]);
			if (L_j > MAXL) continue;
			align_j_flat(labels[j], profiles[j].data(), L_j);
			align_j_old(labels[j], profiles[j].data(), L_j);
			}
		}
	}
