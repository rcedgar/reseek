#include "myutils.h"
#include "flat_bench.h"
#include "alpha.h"
#include "sort.h"
#include "triangle.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "flat_alignx.h"
#include "flat_aligner.h"

static const uint M = 64;

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
	const string &logoddsfnpattern,
	bool set_self_scores)
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
	if (set_self_scores)
		set_selfrev_scores();
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
	float score = fa.m_score;
	//float score2 = calc_aln_weights(fa, DomIdxQ, DomIdxT);//@@TODO

	fa.write_aln(stdout);
	fa.write_aln(g_fLog);
	}

void flat_bench::align_pair_selfrev(FILE *f, uint DomIdx)
	{
	if (f == 0)
		return;

	flat_aligner fa;
	fa.m_ff = &m_ff;
	fa.alloc();

	const uint8_t *profT = m_fp.get_profile(DomIdx);
	const string &label = m_look->get_dom(DomIdx);
	const uint LT = m_fp.get_length(DomIdx);
	fa.cacheT(label, profT, LT);

	uint8_t *profQ = m_fp.get_rev_profile(DomIdx);
	const uint LQ = m_fp.get_length(DomIdx);
	fa.alignQ(label + "_rev", profQ, LQ);
	fprintf(f, "%s\t%.4g\n", label.c_str(), fa.m_score);
	myfree(profQ);
	fa.freemem();
	}

void flat_bench::ThreadBody_All(uint ThreadIdx)
	{
	asserta(!optset_reverse);
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
		if (flat_params::m_rev_w > 0)
			fa.cache_reverseT(labelT, profT, LT);

		// Includes self-score for santify checking and because
		//   triangle*() functions include diagonal
		for (uint DomIdxQ = DomIdxT; DomIdxQ < NQ; ++DomIdxQ)
			{
			const string &labelQ = m_fp.get_label(DomIdxQ);
			const uint8_t *profQ = m_fp.get_profile(DomIdxQ);
			const uint LQ = m_fp.get_length(DomIdxQ);
			fa.alignQ(labelQ, profQ, LQ);
			float Score = fa.m_score;
			if (flat_params::need_reverse())
				fa.align_reverse();
			asserta(!isnan(Score));
			asserta(!isinf(Score));
			uint PairIdx = triangle_ij_to_k(DomIdxT, DomIdxQ, NQ);
			uint progress_count = m_progress_counter++;
#if SHOW_PROGRESS
			if (progress_count%1000 == 0)
				ProgressStep(progress_count, PairCount, "Aligning");
#endif
			const sid_t *distmxQ = 0;
			const sid_t *distmxT = 0;
			float selfQ = FLT_MAX;
			float selfT = FLT_MAX;
			if (flat_params::need_distmx())
				{
				assert(!m_distmxs.empty());
				distmxQ = m_distmxs[DomIdxQ];
				distmxT = m_distmxs[DomIdxT];
				}
			Score += flat_alignx::alignx(
				fa, profQ, profT, distmxQ, distmxT, selfT, selfQ, M);
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
	const float open = -flat_params::m_open;
	const float ext = -flat_params::m_ext;
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
			if (flat_params::m_rev_w > 0)
				fa.cache_reverseT(labelT, profT, LT);
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

//		Score += calc_aln_weights(fa, DomIdxQ, DomIdxT);
		AppendHit(DomIdxT, DomIdxQ, Score);
		}
	}

void flat_bench::Search(const string &how)
	{
	Alloc();

#if SHOW_PROGRESS
	if (how == "all")
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
	else if (how == "dope")
		SB->ThreadBody_Dope(ThreadIdx);
	else
		Die("how=%s", how.c_str());
	}

void flat_bench::SetScalarParams(
	const vector<string> &Names,
	const vector<float> &Values)
	{
	flat_params::set_params(Names, Values);
	}

void flat_bench::ClassifyParams(
	const vector<string> &Names,
	const vector<float> &Values,
	vector<string> &AlphaNames,
	vector<float> &Weights,
	vector<string> &ScalarNames,
	vector<float> &ScalarValues,
	float &selfw,
	float &revw,
	bool &need_distmxs)
	{
	selfw = 0;
	revw = 0;
	need_distmxs = false;

	for (uint i = 0; i < SIZE(Names); ++i)
		{
		const string &Name = Names[i];
		float Value = Values[i];
		if (Name == "open" \
			|| Name == "ext" \
			|| Name == "gap2" \
			|| Name == "selfw" \
			|| Name == "dali" \
			|| Name == "dalix" \
			|| Name == "lddt" \
			|| Name == "entropy" \
			|| Name == "revw")
			{
			ScalarNames.push_back(Name);
			ScalarValues.push_back(Value);
			if (Name == "selfw")
				selfw = Value;
			else if (Name == "revw")
				revw = Value;

			if (Name == "dali" \
				|| Name == "dalix" \
				|| Name == "lddt" \
				|| Name == "entropy")
				need_distmxs = true;
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
	ProgressLog("open=%.3g;", flat_params::m_open);
	ProgressLog("ext=%.3g;", flat_params::m_ext);
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
	float selfw = 0;
	float revw = 0;
	bool need_distmxs;
	flat_bench::ClassifyParams(
		Names, Values, AlphaNames, Weights, ScalarNames, ScalarValues,
		selfw, revw, need_distmxs);
	flat_params::m_rev_w = revw;
	flat_params::m_self_w = selfw;

	if (need_distmxs) asserta(!m_distmxs.empty());

	SetScalarParams(ScalarNames, ScalarValues);

	uint n = SIZE(AlphaNames);
	asserta(SIZE(Weights) == n);
	unordered_map<string, float> NameToWeight;
	for (uint i = 0; i < n; ++i)
		NameToWeight[AlphaNames[i]] = Weights[i];
	ApplyWeightsToLogOdds(NameToWeight);

	if (selfw != 0)
		set_selfrev_scores();
	}

void flat_bench::set_selfrev_scores()
	{
	if (m_self_rev_scores != 0)
		return;
	uint ndom = m_look->get_ndom();
	m_self_rev_scores = myalloc(float, ndom);
	flat_aligner fa;
	fa.m_ff = &m_ff;
	fa.alloc();
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		const string &label = m_look->get_dom(domidx);
		uint L = m_fp.get_length(domidx);
		const uint8_t *prof = m_fp.get_profile(domidx);
		m_self_rev_scores[domidx] =
			fa.get_self_rev_score(label, prof, L);
		}
	fa.freemem();
	}

void flat_bench::set_distmxs(const string &chainfn)
	{
	vector<flat_chain_t *> chains;
	read_flat_chains(chainfn, chains);
	int nchain = int(chains.size());
	const uint ndom = m_look->get_ndom();
	m_distmxs.clear();
	m_distmxs.resize(ndom);
	for (int i = 0; i < nchain; ++i)
		{
		flat_chain_t *chain = chains[i];
		const string &label = chain->m_label;
		uint domidx = m_look->get_domidx(label, true);
		if (domidx == UINT_MAX)
			continue;
		uint L = chain->get_length();
		sid_t *distmx = myalloc(sid_t, L*M);
		chaq::fill_distmx(chain, M, distmx);
		m_distmxs[domidx] = distmx;
		}
	}

void cmd_flat_bench()
	{
	asserta(optset_lookup);
	asserta(optset_fapattern);
	asserta(optset_mxpattern);

	asserta(!optset_spec);
	asserta(!optset_varstr);

	const string &VarStr = g_Arg1;

	flat_bench FB;
	FB.ReadLookup(opt(lookup));
	if (optset_dope)
		FB.ReadDope(opt(dope));
	if (optset_input)
		FB.set_distmxs(opt(input));

	vector<string> param_names;
	vector<float> param_values;
	ParseVarStr(VarStr, param_names, param_values);

	vector<string> feature_names;
	vector<string> scalar_names;
	vector<float> weights;
	vector<float> scalar_values;
	float selfw = 0;
	float revw = 0;
	bool need_distmxs = false;
	flat_bench::ClassifyParams(param_names, param_values,
		feature_names, weights,
		scalar_names, scalar_values,
		selfw, revw, need_distmxs);
	if (need_distmxs)
		asserta(optset_input);

	bool set_self_scores = (selfw != 0);

	FB.load_alphas_and_profiles(
		feature_names, weights, opt(fapattern), opt(mxpattern),
		set_self_scores);
	FB.UpdateParamsFromVarStr(VarStr);
	FB.ProgressLogParams();
	FB.Alloc();
	if (set_self_scores)
		FB.set_selfrev_scores();
	//FB.m_self_rev_weight = selfw;
	//FB.m_rev_weight = revw;
	FB.SetScalarParams(scalar_names, scalar_values);

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
	FB.WriteHits(opt(output), true, opt(triangle));
	}
