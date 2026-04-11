#include "myutils.h"
#include "flat_bench.h"
#include "alpha.h"
#include "sort.h"
#include "triangle.h"
#include "flat_helpers.h"
#include "flat_aligner.h"
#include <unordered_map>
#include <unordered_set>

#define	SHOW_PROGRESS	1

atomic<uint> flat_bench::m_progress_counter;
atomic<uint> flat_bench::m_ncachehits;
atomic<uint> flat_bench::m_ncachemisses;

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values);

void flat_bench::load_alphas_and_profiles(
	const vector<string> &feature_names,
	const vector<float> &weights,
	const string &fafnpattern,
	const string &logoddsfnpattern)
	{
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
	m_fp.read_profiles_from_fastas(fafns, m_look->m_dom2idx);
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
	}

void flat_bench::ThreadBody_All(uint ThreadIdx)
	{
	const uint NQ = SIZE(m_Labels);
	const uint PairCount = triangle_get_K(NQ);
	const uint nfeat = get_nfeat();
	flat_aligner fa;
	fa.m_ff = &m_ff;
	fa.m_open = -m_Open;
	fa.m_ext = -m_Ext;
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
#if SHOW_PROGRESS
			if (progress_count%1000 == 0)
				ProgressStep(progress_count, PairCount, "Aligning");
#endif
			AppendHit(DomIdxT, DomIdxQ, Score);
			}
		}
	}

// Outputs tsv for TS training with CIGAR
// Identical to ThreadBody_All() except for
//			m_lock_tsv_all_vs_all.lock();
//			fa.write_tsv(m_f_tsv_all_vs_all);
//			m_lock_tsv_all_vs_all.unlock();
void flat_bench::ThreadBody_AllVsAll(uint ThreadIdx)
	{
	asserta(m_f_tsv_all_vs_all != 0);
	const uint NQ = SIZE(m_Labels);
	const uint PairCount = triangle_get_K(NQ);
	const uint nfeat = get_nfeat();
	flat_aligner fa;
	fa.m_ff = &m_ff;
	fa.m_open = -m_Open;
	fa.m_ext = -m_Ext;
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

			m_lock_tsv_all_vs_all.lock();
			fa.write_tsv(m_f_tsv_all_vs_all);
			m_lock_tsv_all_vs_all.unlock();

			float Score = fa.m_score;
			asserta(!isnan(Score));
			asserta(!isinf(Score));
			uint PairIdx = triangle_ij_to_k(DomIdxT, DomIdxQ, NQ);
			uint progress_count = m_progress_counter++;
#if SHOW_PROGRESS
			if (progress_count%1000 == 0)
				ProgressStep(progress_count, PairCount, "Aligning");
#endif
			AppendHit(DomIdxT, DomIdxQ, Score);
			}
		}
	}

void flat_bench::ThreadBody_Dope(uint ThreadIdx)
	{
	assert(m_look);
	const uint ndom = m_look->get_ndom();
	const uint nfeat = get_nfeat();
	flat_aligner fa;
	fa.m_ff = &m_ff;
	fa.m_open = -m_Open;
	fa.m_ext = -m_Ext;
	fa.alloc();
	uint CurrentDomIdxT = UINT_MAX;
	for (;;)
		{
		uint dopeidx = m_NextDopeIdx++;
		if (dopeidx >= m_dope_nhit)
			{
			fa.freemem();
			return;
			}
#if SHOW_PROGRESS
		if (dopeidx%1000 == 0)
			ProgressStep(dopeidx, m_dope_nhit, "Aligning");
#endif
		uint k = m_dope_ks[dopeidx];
		uint DomIdxQ, DomIdxT;
		triangle_k_to_ij(k, ndom, DomIdxT, DomIdxQ);

		if (DomIdxT == CurrentDomIdxT)
			++m_ncachehits;
		else
			{
			++m_ncachemisses;
			const string &labelT = m_fp.get_label(DomIdxT);
			const uint8_t *profT = m_fp.get_profile(DomIdxT);
			const uint LT = m_fp.get_length(DomIdxT);
			fa.cacheT(labelT, profT, LT);
			CurrentDomIdxT = DomIdxT;
			}

		const string &labelQ = m_fp.get_label(DomIdxQ);
		const uint8_t *profQ = m_fp.get_profile(DomIdxQ);
		const uint LQ = m_fp.get_length(DomIdxQ);
		fa.alignQ(labelQ, profQ, LQ);
		float Score = fa.m_score;
		asserta(!isnan(Score));
		asserta(!isinf(Score));
		uint PairIdx = triangle_ij_to_k(DomIdxT, DomIdxQ, ndom);
		uint progress_count = m_progress_counter++;
		AppendHit(DomIdxT, DomIdxQ, Score);
		}
	}

void flat_bench::Search(const string &how)
	{
	Alloc();

#if SHOW_PROGRESS
	if (how == "all" || how == "allvsall")
		{
		const uint NQ = SIZE(m_Labels);
		const uint PairCount = triangle_get_K(NQ);
		ProgressStep(0, PairCount, "Aligning");
		}
	else if (how == "dope")
		{
		ProgressStep(0, m_dope_nhit, "Aligning");
		}
#endif

	m_ThreadCount = GetRequestedThreadCount();
	m_NextQueryIdx = 0;
	m_NextDopeIdx = 0;
	m_progress_counter = 0;
	m_ncachehits = 0;
	m_ncachemisses = 0;
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadBody, this, ThreadIndex, how);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];

#if SHOW_PROGRESS
	if (how == "all")
		{
		const uint NQ = SIZE(m_Labels);
		const uint PairCount = triangle_get_K(NQ);
		ProgressStep(PairCount-1, PairCount, "Aligning");
		}
	else if (how == "dope")
		{
		ProgressStep(m_dope_nhit-1, m_dope_nhit, "Aligning");
		uint hits = m_ncachehits;
		uint misses = m_ncachemisses;
		ProgressLog("Cache misses %u, hits %u (%.1f%%)\n",
			misses, hits, GetPct(hits, hits+misses));
		}
#endif
	}

void flat_bench::StaticThreadBody(flat_bench *SB,
	uint ThreadIdx, const string &how)
	{
	if (how == "all")
		SB->ThreadBody_All(ThreadIdx);
	if (how == "allvsall")
		SB->ThreadBody_AllVsAll(ThreadIdx);
	else if (how == "dope")
		SB->ThreadBody_Dope(ThreadIdx);
	else
		Die("how=%s", how.c_str());
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
	asserta(optset_lookup);
	asserta(optset_dope);
	asserta(optset_fapattern);
	asserta(optset_mxpattern);

	asserta(!optset_spec);
	asserta(!optset_varstr);

	const string &VarStr = g_Arg1;

	flat_bench FB;
	FB.ReadLookup(opt(lookup));
	if (optset_dope)
		FB.ReadDope(opt(dope));

	vector<string> param_names;
	vector<float> param_values;
	ParseVarStr(VarStr, param_names, param_values);

	vector<string> feature_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	flat_bench::ClassifyParams(param_names, param_values,
		feature_names, weights,
		scalar_names, scalar_values);

	FB.load_alphas_and_profiles(
		feature_names, weights, opt(fapattern), opt(mxpattern));
	FB.UpdateParamsFromVarStr(VarStr);
	FB.ProgressLogParams();
	FB.Alloc();

	if (optset_label1)
		{
		asserta(optset_label2);
		FB.align_pair(opt(label1), opt(label2));
		return;
		}

	if (optset_dope)
		FB.Search("dope");
	else
		FB.Search("all");
	FB.SetScoreOrder();
	FB.Bench();
	FB.WriteHits(opt(output), true);
	}
