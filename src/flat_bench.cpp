#include "myutils.h"
#include "flat_bench.h"
#include "alpha.h"
#include "sort.h"
#include "triangle.h"
#include "flat_helpers.h"
#include "flat_aligner.h"
#include <unordered_map>
#include <unordered_set>

atomic<uint> flat_bench::m_progress_counter;

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values);

void flat_bench::load_alphas_and_profiles(
	const string &VarStr,
	const string &fafnpattern,
	const string &logoddsfnpattern)
	{
	vector<string> param_names;
	vector<float> param_values;
	ParseVarStr(VarStr, param_names, param_values);

	vector<string> feature_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	ClassifyParams(param_names, param_values,
		feature_names, weights,
		scalar_names, scalar_values);

	const uint nfeat = uint(feature_names.size());
	asserta(weights.size() == nfeat);

	vector<string> fafns(nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		make_fn_pattern(
			fafnpattern,
			feature_names[fi],
			fafns[fi]);

	m_ff.init(feature_names);
	m_fp.m_ff = &m_ff;
	m_fp.read_profiles_from_fastas(fafns);
	m_ff.read_logoddsvec_pattern(logoddsfnpattern);
	m_ff.set_feature_block_offsets();
	m_ff.set_symbolsvec();

	m_Labels = m_fp.m_labels;
	m_SeqCount = uint(m_Labels.size());
	}

void flat_bench::align_pair(
	const string &labelQ, const string &labelT)
	{
	flat_aligner fa;
	fa.m_ff = &m_ff;
	fa.alloc();

	uint DomIdxQ = UINT_MAX;
	uint DomIdxT = UINT_MAX;
	for (uint i = 0; i < uint(m_Labels.size()); ++i)
		{
		if (m_Labels[i] == labelQ)
			DomIdxQ = i;
		if (m_Labels[i] == labelT)
			DomIdxT = i;
		}
	asserta(DomIdxQ != UINT_MAX);
	asserta(DomIdxT != UINT_MAX);

	m_fp.profile_to_fasta(g_fLog, DomIdxQ);
	m_fp.profile_to_fasta(g_fLog, DomIdxT);

	const uint8_t *profT = m_fp.get_profile(DomIdxT);
	const uint LT = m_fp.get_length(DomIdxT);
	fa.cacheT(labelT, profT, LT);

	const uint8_t *profQ = m_fp.get_profile(DomIdxQ);
	const uint LQ = m_fp.get_length(DomIdxQ);
	fa.alignQ(labelQ, profQ, LQ);

	fa.write_aln(stdout);
	fa.write_aln(g_fLog);

	bool tp = IsTP(DomIdxQ, DomIdxT);
	ProgressLog("%s\n", tp ? "TRUE POSITIVE" : "FALSE POSITIVE");
	}

void flat_bench::ThreadBody_All(uint ThreadIdx)
	{
	const uint NQ = SIZE(m_Labels);
	const uint PairCount = triangle_get_K(NQ);
	const uint nfeat = get_nfeat();
	flat_aligner fa;
	fa.m_ff = &m_ff;
	fa.alloc();
	for (;;)
		{
		uint DomIdxT = m_NextQueryIdx++;
		if (DomIdxT >= NQ)
			return;

		const string &labelT = m_fp.get_label(DomIdxT);
		const uint8_t *profT = m_fp.get_profile(DomIdxT);
		const uint LT = m_fp.get_length(DomIdxT);
		fa.cacheT(labelT, profT, LT);

		// Includes self-score for santify checking and because
		//   triangle*() functions include diagonal
		for (uint DomIdxQ = DomIdxT; DomIdxQ < NQ; ++DomIdxQ)
			{
			const string &labelQ = m_fp.get_label(DomIdxQ);
			const uint8_t *profQ = m_fp.get_profile(DomIdxQ);
			const uint LQ = m_fp.get_length(DomIdxQ);
			fa.alignQ(labelQ, profQ, LQ);
			float Score = fa.m_score;
			asserta(!isnan(Score));
			asserta(!isinf(Score));
			uint PairIdx = triangle_ij_to_k(DomIdxT, DomIdxQ, NQ);
			uint progress_count = m_progress_counter++;
			if (progress_count%1000 == 0)
				ProgressStep(progress_count, PairCount, "Aligning");
			AppendHit(DomIdxT, DomIdxQ, Score);
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
	const uint PairCount = triangle_get_K(NQ);
	ProgressStep(PairCount-1, PairCount, "Aligning");
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
	m_ff.apply_weights(NameToWeight);
	}

void flat_bench::ProgressLogParams() const
	{
	ProgressLog("open=%.3g;", m_Open);
	ProgressLog("ext=%.3g;", m_Ext);
	uint nfeat = m_ff.get_nfeat();
	for (uint fi = 0; fi < nfeat; ++fi)
		ProgressLog("%s=%.3g;",
			m_ff.m_feature_names[fi].c_str(),
			m_ff.m_weights[fi]);
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

void cmd_flat_bench()
	{
	asserta(!optset_lookup);
	asserta(!optset_spec);
	asserta(!optset_varstr);
	asserta(optset_fapattern);
	asserta(optset_mxpattern);

	const string &VarStr = g_Arg1;

	flat_bench FB;
	FB.load_alphas_and_profiles(
		VarStr, opt(fapattern), opt(mxpattern));
	FB.SetLookupFromLabels();
	FB.UpdateParamsFromVarStr(VarStr);
	FB.ProgressLogParams();
	FB.Alloc();

	if (optset_label1)
		{
		asserta(optset_label2);
		FB.align_pair(opt(label1), opt(label2));
		return;
		}

	FB.Search_All();
	FB.SetScoreOrder();
	FB.Bench_All();
	FB.WriteHits(opt(output), true);
	}
