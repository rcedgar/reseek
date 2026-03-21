#include "myutils.h"
#include "dss.h"
#include "chaq.h"
#include "flat_base.h"
#include "flat_chain.h"
#include "pdbchain.h"
#include "flat_distmx.h"
#include "pdbfilescanner.h"
#include "flat_chain_reader.h"
#include "flat_helpers.h"
#include "flat_features.h"
#include "flat_profiles.h"
#include "flat_aligner.h"
#include "getticks.h"

static const uint s_maxL = 4000;

static uint s_nfeat;
static uint *s_alpha_sizes;
static uint *s_feature_block_offsets;
static float *s_pssm_i;
static float **s_weighted_logoddsmxvec;
static const uint8_t *s_prof_i;
static uint s_L_i;
static string s_label_i;

static const float **__restrict s_scratch_pssms;
static float *s_scratch_rows;
static uint8_t *s_TB;

static float s_open = -3;
static float s_ext = -1;

static void cache_i(const string &label, const uint8_t *prof_i, uint L_i)
	{
	s_label_i = label;
	s_prof_i = prof_i;
	s_L_i = L_i;
	fill_flat_pssm(prof_i, L_i, s_nfeat, s_alpha_sizes,
		s_feature_block_offsets, s_weighted_logoddsmxvec, s_pssm_i);
	}

static float align_j(
	const string &label,
	const uint8_t *prof_j,
	uint L_j,
	uint &Loi,
	uint &Loj,
	char *path_buffer,
	uint &ncol)
	{
	float score = sw_flat_pssm(
		s_scratch_rows, s_TB, s_scratch_pssms,
		prof_j, L_j,
		s_pssm_i, s_L_i, s_feature_block_offsets,
		s_nfeat, s_open, s_ext,
		Loi, Loj, path_buffer, ncol);
	return score;
	}

static float align_j(const string &label, const uint8_t *prof_j, uint L_j)
	{
	uint Loi, Loj;
	char *path_buffer = myalloc(char, 2*s_maxL);
	uint ncol;
	float score = align_j(label, prof_j, L_j,
		Loi, Loj, path_buffer, ncol);
	return score;
	}

void cmd_flat_align_pairs_spec()
	{
	const string &specfn = g_Arg1;
	vector<string> labels;
	vector<vector<uint8_t> > profiles;
	vector<string> feature_names;
	vector<uint> alpha_sizes;
	vector<vector<float> > logoddsmxvec;
	read_profiles_and_logoddsvec(
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
	const uint32_t sum_alpha_sizes =
		get_flat_pssm_feature_block_offsets(s_nfeat,
			s_alpha_sizes, s_feature_block_offsets);

	s_pssm_i = myalloc(float, s_maxL * sum_alpha_sizes);

	s_scratch_rows = myalloc(float, 2*s_maxL + 2);
	s_scratch_pssms = myalloc(const float *, s_nfeat);
	s_TB = myalloc(uint8_t, s_maxL*s_maxL);

	uint npairs = nprof*nprof;

	ProgressLog("%10u  features\n", s_nfeat);
	ProgressLog("%10u  profiles\n", nprof);
	ProgressLog("%10u  pairs\n", npairs);

	uint counter = 0;
	TICKS t1 = GetClockTicks();
	for (uint i = 0; i < nprof; ++i)
		{
		uint L_i = SIZE(profiles[i]);
		asserta(L_i%s_nfeat == 0);
		L_i /= s_nfeat;
		if (L_i > s_maxL) continue;
		cache_i(labels[i], profiles[i].data(), L_i);

		for (uint j = 0; j < nprof; ++j)
			{
			//ProgressStep(counter++, npairs, "Aligning");
			uint L_j = SIZE(profiles[j]);
			asserta(L_j%s_nfeat == 0);
			L_j /= s_nfeat;
			if (L_j > s_maxL) continue;
			align_j(labels[j], profiles[j].data(), L_j);
			}
		}
	TICKS t2 = GetClockTicks();
	double t = double(t2 - t1);
	ProgressLog("%.3g ticks\n", t);
	}

void cmd_flat_align_pairs_faprof()
	{
	asserta(optset_mxpattern);
	const string &faproffn = g_Arg1;

	flat_profiles fp;
	vector<string> feature_names;
	vector<uint> alpha_sizes;
	fp.read_profiles_faprof(faproffn, feature_names);
	fp.m_ff = new flat_features;
	fp.m_ff->init(feature_names);
	fp.m_ff->read_logoddsvec_pattern(opt(mxpattern));
	fp.m_ff->apply_unit_weights();
	fp.check_profiles();

	flat_aligner fa;
	fa.m_ff = fp.m_ff;
	fa.alloc();

	const uint nprof = fp.get_nprof();
	const uint npairs = nprof*nprof;

	ProgressLog("%10u  features\n", s_nfeat);
	ProgressLog("%10u  profiles\n", nprof);
	ProgressLog("%10u  pairs\n", npairs);

	for (uint i = 0; i < nprof; ++i)
		{
		uint L_i = fp.get_length(i);
		if (L_i > s_maxL) continue;
		const string &label_i = fp.get_label(i);
		const uint8_t *prof_i = fp.get_profile(i);
		fa.cacheT(label_i, prof_i, L_i);

		for (uint j = 0; j < nprof; ++j)
			{
			uint L_j = fp.get_length(j);
			if (L_j > s_maxL) continue;
			const string &label_j = fp.get_label(j);
			const uint8_t *prof_j = fp.get_profile(j);
			fa.alignQ(label_j, prof_j, L_j);
			fa.write_aln(g_fLog);
			}
		}
	}
