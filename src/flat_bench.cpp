#include "myutils.h"
#include "flat_bench.h"
#include "alpha.h"
#include "sort.h"
#include "triangle.h"
#include "flat_helpers.h"
#include <unordered_map>
#include <unordered_set>

atomic<uint> flat_bench::m_progress_counter;

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values);

void flat_bench::load_alphas_and_profiles(const string &SpecFN)
	{
	read_profiles_and_logoddsvec(
		SpecFN,
		m_AlphaNames,
		m_alpha_sizes,
		m_Labels,
		m_profiles,
		m_raw_logoddsmxvec);


	uint nfeat = SIZE(m_AlphaNames);
	asserta(SIZE(m_alpha_sizes) == nfeat);
	m_weighted_logoddsmxvec.clear();
	m_weighted_logoddsmxvec.resize(nfeat);

// Initialize to uniform weights
	const float w = 1.0f/nfeat;;
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = m_alpha_sizes[fi];
		m_weighted_logoddsmxvec[fi].reserve(AS);
		asserta(SIZE(m_raw_logoddsmxvec[fi]) == AS*AS);
		for (uint k = 0; k < AS*AS; ++k)
			m_weighted_logoddsmxvec[fi].push_back(
				w*m_raw_logoddsmxvec[fi][k]);
		}

	uint NQ = SIZE(m_Labels);
	asserta(SIZE(m_profiles) == NQ);

	m_MaxL = 0;
	for (uint i = 0; i < NQ; ++i)
		{
		uint profile_length = SIZE(m_profiles[i]);
		asserta(profile_length%nfeat == 0);
		uint L = profile_length/nfeat;
		m_MaxL = max(m_MaxL, L);
		}
	ProgressLog("Max length %u\n", m_MaxL);
	}

void flat_bench::ThreadBody_All(uint ThreadIdx)
	{
	const uint NQ = SIZE(m_Labels);
	const uint PairCount = triangle_get_K(NQ);
	const uint nfeat = SIZE(m_AlphaNames);
	asserta(SIZE(m_weighted_logoddsmxvec) == nfeat);
	float **weighted_logoddsmxvec = myalloc(float *, nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = m_alpha_sizes[fi];
		asserta(AS >= 2 && AS < 256);
		weighted_logoddsmxvec[fi] = m_weighted_logoddsmxvec[fi].data();
		}
	uint32_t *feature_block_offsets = myalloc(uint32_t, nfeat);
	const uint32_t sum_alpha_sizes =
		get_flat_pssm_feature_block_offsets(nfeat,
			m_alpha_sizes.data(), feature_block_offsets);

	float *pssmQ = myalloc(float, m_MaxL*sum_alpha_sizes);

	float *scratch_rows = myalloc(float, 2*m_MaxL + 2);
	const float **scratch_pssms = myalloc(const float *, nfeat);
	uint8_t *TB = myalloc(uint8_t, m_MaxL*m_MaxL);
	uint Loi, Loj;
	string Path;
	for (;;)
		{
		uint DomIdxQ = m_NextQueryIdx++;
		if (DomIdxQ >= NQ)
			return;

		/////////////////////////////////////////////////////////
		// Cache query
		/////////////////////////////////////////////////////////
		const vector<uint8_t> &profvecQ = m_profiles[DomIdxQ];
		uint profile_length = SIZE(profvecQ);
		asserta(profile_length%nfeat == 0);
		uint LQ = profile_length/nfeat;
		assert(LQ <= m_MaxL);
		if (SIZE(profvecQ) != LQ*nfeat)
			{
			ProgressLog("DomIdxQ       %u\n", DomIdxQ);
			ProgressLog("nfeat      %u\n", nfeat);
			ProgressLog("Label %s\n", m_Labels[DomIdxQ].c_str());
			ProgressLog("LQ         %u\n", LQ);
			ProgressLog("proflen    %u\n", SIZE(profvecQ));
			ProgressLog("LQ*nfeat   %u\n", LQ*nfeat);
			Die("SIZE(profvecQ) != LQ*nfeat");
			}
		const uint8_t *profQ = profvecQ.data();
		fill_flat_pssm(profQ, LQ, nfeat, m_alpha_sizes.data(),
			feature_block_offsets, weighted_logoddsmxvec, pssmQ);
		/////////////////////////////////////////////////////////

		// Includes self-score for santify checking and because
		//   triangle*() functions include diagonal
		for (uint DomIdxT = DomIdxQ; DomIdxT < NQ; ++DomIdxT)
			{
			const vector<uint8_t> &profvecT = m_profiles[DomIdxT];
			uint profile_lengthT = SIZE(profvecT);
			assert(profile_lengthT%nfeat == 0);
			uint LT = profile_lengthT/nfeat;
			assert(LT <= m_MaxL);
			const uint8_t *profT = profvecT.data();

			float Score = sw_flat_pssm(
				scratch_rows, TB, scratch_pssms,
				profT, LT,
				pssmQ, LQ,
				feature_block_offsets, nfeat,
				-m_Open, -m_Ext,
				Loi, Loj, Path);

			asserta(!isnan(Score));
			asserta(!isinf(Score));
			uint PairIdx = triangle_ij_to_k(DomIdxQ, DomIdxT, NQ);
			uint progress_count = m_progress_counter++;
			if (progress_count%1000 == 0)
				ProgressStep(progress_count, PairCount, "Aligning");

			AppendHit(DomIdxQ, DomIdxT, Score);
			}
		}
	}

void flat_bench::Search_All()
	{
	m_ThreadCount = GetRequestedThreadCount();
	m_NextQueryIdx = 0;
	m_progress_counter = 0;
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadBody_All, this, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	const uint NQ = SIZE(m_Labels);
	ProgressStep(NQ-1, NQ, "Aligning");
	}

void flat_bench::StaticThreadBody_All(flat_bench *SB, uint ThreadIdx)
	{
	SB->ThreadBody_All(ThreadIdx);
	}

void flat_bench::Bench_All(const string &Msg)
	{
	SetScoreOrder();
	FastBench::Bench(Msg);
	}

void flat_bench::SetScalarParams(
	const vector<string> &Names,
	const vector<float> &Values)
	{
	m_Open = FLT_MAX;
	m_Ext = FLT_MAX;
	for (uint i = 0; i < SIZE(Names); ++i)
		{
		const string &Name = Names[i];
		float Value = Values[i];
		if (Name == "open")
			m_Open = Value;
		else if (Name == "ext")
			m_Ext = Value;
		else if (Name == "gap2")
			{
			m_Open = Value;
			m_Ext = Value/10;
			}
		else
			Die("flat_bench::SetScalarParams() Name=%s", Name.c_str());
		}
	asserta(m_Open != FLT_MAX && m_Ext != FLT_MAX);
	}

void flat_bench::ClassifyParams(
	const vector<string> &Names,
	const vector<float> &Values,
	vector<string> &AlphaNames,
	vector<float> &Weights,
	vector<string> &ScalarNames,
	vector<float> &ScalarValues)
	{
	for (uint i = 0; i < SIZE(Names); ++i)
		{
		const string &Name = Names[i];
		float Value = Values[i];
		if (Name == "open" || Name == "ext" || Name == "gap2")
			{
			ScalarNames.push_back(Name);
			ScalarValues.push_back(Value);
			}
		else
			{
			AlphaNames.push_back(Name);
			Weights.push_back(Value);
			}
		}
	}

void flat_bench::ApplyWeightsToLogOdds(
	const unordered_map<string, float> &NameToWeight)
	{
	uint nalpha = SIZE(m_AlphaNames);
	m_Weights.clear();
	m_Weights.resize(nalpha, 0);
	asserta(SIZE(NameToWeight) == nalpha);
	unordered_map<string, uint> NameToIdx;
	for (uint idx = 0; idx < nalpha; ++idx)
		NameToIdx[m_AlphaNames[idx]] = idx;

	for (unordered_map<string, float>::const_iterator iter = NameToWeight.begin();
		iter != NameToWeight.end(); ++iter)
		{
		const string &Name = iter->first;
		float Weight = iter->second;
		unordered_map<string, uint>::const_iterator iter2 = NameToIdx.find(Name);
		asserta(iter2 != NameToIdx.end());
		uint idx = iter2->second;
		m_Weights[idx] = Weight;

		uint AS = m_alpha_sizes[idx];
		for (uint code = 0; code < AS; ++code)
			m_weighted_logoddsmxvec[idx][code] =
				m_raw_logoddsmxvec[idx][code]*Weight;
		}
	}

void flat_bench::ProgressLogParams() const
	{
	ProgressLog("open=%.3g;", m_Open);
	ProgressLog("ext=%.3g;", m_Ext);
	uint nfeat = SIZE(m_AlphaNames);
	for (uint fi = 0; fi < nfeat; ++fi)
		ProgressLog("%s=%.3g;",
			m_AlphaNames[fi].c_str(), m_Weights[fi]);
	ProgressLog("\n");
	}

void flat_bench::UpdateParamsFromVarStr(const string &VarStr)
	{
	vector<string> Names;
	vector<float> Values;
	ParseVarStr(VarStr, Names, Values);

	vector<string> AlphaNames;
	vector<float> Weights;
	vector<string> ScalarNames;
	vector<float> ScalarValues;
	flat_bench::ClassifyParams(
		Names, Values, AlphaNames, Weights, ScalarNames, ScalarValues);

	SetScalarParams(ScalarNames, ScalarValues);

	uint n = SIZE(AlphaNames);
	asserta(SIZE(Weights) == n);
	unordered_map<string, float> NameToWeight;
	for (uint i = 0; i < n; ++i)
		NameToWeight[AlphaNames[i]] = Weights[i];
	ApplyWeightsToLogOdds(NameToWeight);
	}

void flat_bench::AlignPair(const string &LabelQ, const string &LabelT)
	{
	uint profile_idxQ = UINT_MAX;
	uint profile_idxT = UINT_MAX;
	const uint NQ = SIZE(m_Labels);
	for (uint i = 0; i < NQ; ++i)
		{
		if (m_Labels[i] == LabelQ)
			profile_idxQ = i;
		if (m_Labels[i] == LabelT)
			profile_idxT = i;
		}
	asserta(profile_idxQ < NQ && profile_idxT < NQ);
	const uint nfeat = SIZE(m_AlphaNames);
	asserta(SIZE(m_weighted_logoddsmxvec) == nfeat);
	float **weighted_logoddsmxvec = myalloc(float *, nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		{
		uint AS = m_alpha_sizes[fi];
		asserta(AS >= 2 && AS < 256);
		weighted_logoddsmxvec[fi] = m_weighted_logoddsmxvec[fi].data();
		}
	uint32_t *feature_block_offsets = myalloc(uint32_t, nfeat);
	const uint32_t sum_alpha_sizes =
		get_flat_pssm_feature_block_offsets(nfeat,
			m_alpha_sizes.data(), feature_block_offsets);

	float *pssmQ = myalloc(float, m_MaxL*sum_alpha_sizes);

	float *scratch_rows = myalloc(float, 2*m_MaxL + 2);
	const float **scratch_pssms = myalloc(const float *, nfeat);
	uint8_t *TB = myalloc(uint8_t, m_MaxL*m_MaxL);
	uint Loi, Loj;
	string Path;

	const vector<uint8_t> &profvecQ = m_profiles[profile_idxQ];
	uint profile_length = SIZE(profvecQ);
	asserta(profile_length%nfeat == 0);
	uint LQ = profile_length/nfeat;
	assert(LQ <= m_MaxL);
	asserta(SIZE(profvecQ) == LQ*nfeat);

	const uint8_t *profQ = profvecQ.data();
	fill_flat_pssm(profQ, LQ, nfeat, m_alpha_sizes.data(),
		feature_block_offsets, weighted_logoddsmxvec, pssmQ);

	const vector<uint8_t> &profvecT = m_profiles[profile_idxT];
	uint profile_lengthT = SIZE(profvecT);
	assert(profile_lengthT%nfeat == 0);
	uint LT = profile_lengthT/nfeat;
	assert(LT <= m_MaxL);
	const uint8_t *profT = profvecT.data();

	float Score = sw_flat_pssm(
		scratch_rows, TB, scratch_pssms,
		profT, LT,
		pssmQ, LQ,
		feature_block_offsets, nfeat,
		-m_Open, -m_Ext,
		Loi, Loj, Path);

	Log("score %.3g  %s\n", Score, Path.c_str());
	}

void cmd_flat_bench_align_pair()
	{
	asserta(!optset_lookup);
	asserta(!optset_spec);
	asserta(optset_varstr);
	asserta(optset_label1);
	asserta(optset_label2);

	const string &SpecFN = g_Arg1;
	const string &VarStr = opt(varstr);

	vector<string> Names;
	vector<float> Values;
	ParseVarStr(VarStr, Names, Values);

	flat_bench FB;
	FB.load_alphas_and_profiles(SpecFN);
	FB.SetLookupFromLabels();
	FB.UpdateParamsFromVarStr(VarStr);
	FB.ProgressLogParams();
	FB.AlignPair(opt(label1), opt(label2));
	}

void cmd_flat_bench()
	{
	asserta(!optset_lookup);
	asserta(!optset_spec);
	asserta(optset_varstr);

	const string &SpecFN = g_Arg1;
	const string &VarStr = opt(varstr);

	vector<string> Names;
	vector<float> Values;
	ParseVarStr(VarStr, Names, Values);

	flat_bench FB;
	FB.load_alphas_and_profiles(SpecFN);
	FB.SetLookupFromLabels();
	FB.UpdateParamsFromVarStr(VarStr);
	FB.ProgressLogParams();
	FB.Alloc();
	FB.Search_All();
	FB.SetScoreOrder();
	FB.Bench_All();
	FB.WriteHits(opt(output), true);
	}
