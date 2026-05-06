#include "myutils.h"
#include "flat_bench.h"
#include "alpha.h"
#include "sort.h"
#include "triangle.h"
#include "flat_params.h"
#include "flat_helpers.h"
#include "flat_alignx.h"
#include "flat_aligner.h"
#include "thread_affinity.h"

#define	SHOW_PROGRESS	1

static FILE *s_f3;	// structure feature scores

atomic<uint> flat_bench::m_aligned_pair_count;
atomic<uint> flat_bench::m_ncachehits;
atomic<uint> flat_bench::m_ncachemisses;
atomic<bool> flat_bench::m_max_secs_exceeded;

static FILE *s_ftsv;
static mutex s_ftsv_lock;

void ParseVarStr(
	const string &VarStr,
	vector<string> &Names,
	vector<float> &Values);

void flat_bench::StaticThreadBody_MaxSecs(uint MaxSecs)
	{
	asserta(MaxSecs > 0);
	m_max_secs_exceeded = false;
	if (MaxSecs == UINT_MAX)
		return;
	std::this_thread::sleep_for(std::chrono::seconds(MaxSecs));
	m_max_secs_exceeded = true;
	}

void flat_bench::load_profiles(const string &fafnpattern)
	{
	const uint nfeat = flat_features::get_nfeat();

	vector<string> fafns(nfeat);
	for (uint fi = 0; fi < nfeat; ++fi)
		make_fn_pattern(
			fafnpattern,
			flat_features::m_feature_names[fi],
			fafns[fi]);

	m_fp.read_profiles_from_fastas(fafns, m_look->m_dom2idx);

	m_Labels = m_fp.m_labels;
	m_SeqCount = uint(m_Labels.size());
	}

void flat_bench::align_pair(flat_aligner &fa,
	const string &labelQ, const string &labelT)
	{
	fa.alloc();

	string labQ, labT;
	trunc_label(labelQ, labQ);
	trunc_label(labelT, labT);
	uint DomIdxQ = UINT_MAX;
	uint DomIdxT = UINT_MAX;
	for (uint i = 0; i < uint(m_Labels.size()); ++i)
		{
		if (m_Labels[i] == labQ)
			DomIdxQ = i;
		if (m_Labels[i] == labT)
			DomIdxT = i;
		}
	asserta(DomIdxQ != UINT_MAX);
	asserta(DomIdxT != UINT_MAX);

	//m_fp.profile_to_fasta(g_fLog, DomIdxQ);
	//m_fp.profile_to_fasta(g_fLog, DomIdxT);

	const uint8_t *profT = m_fp.get_profile(DomIdxT);
	const uint LT = m_fp.get_length(DomIdxT);
	fa.cacheT(labelT, profT, LT);

	const uint8_t *profQ = m_fp.get_profile(DomIdxQ);
	const uint LQ = m_fp.get_length(DomIdxQ);
	fa.alignQ(labelQ, profQ, LQ);
	float score = fa.m_score;
	asserta(!flat_params::need_alignx());
	fa.write_aln(g_fLog);
	if (optset_output)
		{
		FILE *f = CreateStdioFile(opt(output));
		fprintf(f, "%.3g\t%s\t%s\n",
			score, labelQ.c_str(), labelT.c_str());
		CloseStdioFile(f);
		}
	}

void flat_bench::align_pair_selfrev(FILE *f, uint DomIdx)
	{
	if (f == 0)
		return;

	flat_aligner fa;
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

void flat_bench::doT(flat_aligner &fa, uint domidxT)
	{
	const string &labelT = m_fp.get_label(domidxT);
	const uint8_t *profT = m_fp.get_profile(domidxT);
	const uint LT = m_fp.get_length(domidxT);
	fa.cacheT(labelT, profT, LT);
	if (flat_params::m_rev_w > 0)
		fa.cache_reverseT(labelT, profT, LT);
	}

void flat_bench::doQ(flat_aligner &fa, uint domidxQ, uint domidxT)
	{
	if (!in_dope(domidxQ, domidxT))
		{
		AppendHit(domidxQ, domidxT, get_missing_score());
		return;
		}

	const string &labelQ = m_fp.get_label(domidxQ);
	const uint8_t *profQ = m_fp.get_profile(domidxQ);
	const uint8_t *profT = m_fp.get_profile(domidxT);
	const uint LQ = m_fp.get_length(domidxQ);
	fa.alignQ(labelQ, profQ, LQ);
	float dpscore = fa.m_score;
	float Score = fa.m_score;
	asserta(!isnan(Score));
	asserta(!isinf(Score));

	const sid_t *distmxQ = 0;
	const sid_t *distmxT = 0;
	float selfQ = FLT_MAX;
	float selfT = FLT_MAX;
	if (flat_params::need_distmx())
		{
		assert(!m_distmxs.empty());
		distmxQ = m_distmxs[domidxQ];
		distmxT = m_distmxs[domidxT];
		}
	if (flat_params::need_self())
		{
		assert(m_self_rev_scores != 0);
		selfQ = m_self_rev_scores[domidxQ];
		selfT = m_self_rev_scores[domidxT];
		}
	if (s_f3 != 0)
		{
		static mutex s_lock;
		asserta(distmxQ && distmxT);
		float lddt = flat_getlddt_muscle_some_floats4(
			fa, distmxQ, distmxT);
		float dali = 
			flat_get_dali3(fa, distmxQ, distmxT);
		uint n = uint(fa.m_ncol);
		float *colscores = myalloc(float, n);
		float dalix = 
			flat_get_dalix3(fa, distmxQ, distmxT, colscores);
		myfree(colscores);
		float LQ = (float) fa.m_LQ;
		float LT = (float) fa.m_LT;
		float L = (LQ + LT)/2.0f + 50;

		uint ncol = uint(fa.m_ncol);
		float Lfactor = ncol/L;
		float lddtx = lddt*Lfactor;

		s_lock.lock();
		fprintf(s_f3, "%s\t%s\t%.3g\t%.3g\t%.3g\n",
			fa.m_labelQ.c_str(), fa.m_labelT.c_str(), lddtx, dali, dalix);
		s_lock.unlock();
		}
	Score = flat_alignx::alignx(
		fa, profQ, profT, distmxQ, distmxT, selfT, selfQ);
	if (flat_params::need_reverse())
		{
		asserta(flat_params::m_rev_w > 0);
		fa.align_reverse();
		if (flat_params::m_rev_w > 0)
			{
			asserta(fa.m_reverse_score_set);
			Score -= flat_params::m_rev_w*fa.m_reverse_score;
			}
		}
	asserta(!isnan(Score));
	asserta(!isinf(Score));

	AppendHit(domidxT, domidxQ, Score);

	if (!s_ftsv) return;

	string qacc, tacc;
	string line;
	trunc_label(fa.m_labelQ, qacc);
	trunc_label(fa.m_labelT, tacc);
	float lddt = flat_getlddt_old(fa, distmxQ, distmxT);
	uint l2 = (LQ + fa.m_LT)/2;
          
	// qacc+tacc+dpscore+selfrevq+selfrevt+selfrev+lddt+l2+newts
	Ps(line, "%s", qacc.c_str());
	Psa(line, "\t%s", tacc.c_str());
	Psa(line, "\t%.3g", dpscore);
	Psa(line, "\t%.3g", selfT);
	Psa(line, "\t%.3g", selfQ);
	Psa(line, "\t%.3g", (selfQ + selfT)/2);
	Psa(line, "\t%.3g", lddt);
	Psa(line, "\t%u", l2);
	Psa(line, "\t%.3g", Score);
	line += '\n';

	s_ftsv_lock.lock();
	fputs(line.c_str(), s_ftsv);
	s_ftsv_lock.unlock();
	}

void flat_bench::Launch(
	uint ThreadCount, bool PinThreads, bool UseDope, uint MaxSecs)
	{
	asserta(MaxSecs != 0);
	m_NextQueryIdx = 0;
	m_NextDopeIdx = 0;
	m_aligned_pair_count = 0;
	m_ncachehits = 0;
	m_ncachemisses = 0;

	thread *max_secs_thread = 0;
	if (MaxSecs != UINT_MAX)
		max_secs_thread = new thread(StaticThreadBody_MaxSecs, MaxSecs);

	thread_affinity ta;
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(StaticThreadBody, this, ThreadIndex, UseDope);
		if (PinThreads)
			ta.pinThread(*t, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		delete ts[ThreadIndex];
	if (max_secs_thread != 0)
		{
		max_secs_thread->join();
		delete max_secs_thread;
		}
	}


void flat_bench::ThreadBody_All(uint ThreadIdx)
	{
	m_aligned_pair_count = 0;
	const uint NQ = SIZE(m_Labels);
	const uint PairCount = triangle_get_K(NQ);
	const uint nfeat = flat_features::get_nfeat();
	flat_aligner fa;
	fa.alloc();
	for (;;)
		{
		uint DomIdxT = m_NextQueryIdx++;
		if (DomIdxT >= NQ)
			{
			fa.freemem();
			return;
			}

		doT(fa, DomIdxT);
		// Includes self-score for santify checking and because
		//   triangle*() functions include diagonal
		for (uint DomIdxQ = DomIdxT; DomIdxQ < NQ; ++DomIdxQ)
			{
			if (m_max_secs_exceeded)
				{
				fa.freemem();
				return;
				}
			uint PairIdx = triangle_ij_to_k(DomIdxT, DomIdxQ, NQ);
			doQ(fa, DomIdxQ, DomIdxT);
			uint progress_count = m_aligned_pair_count++;
#if SHOW_PROGRESS
			if (ThreadIdx == 0 && progress_count%1000 == 0)
				ProgressStep(progress_count, PairCount, "Aligning");
#endif
			}
		}
	fa.freemem();
	}

void flat_bench::ThreadBody_Dope(uint ThreadIdx)
	{
	assert(m_look);
	m_aligned_pair_count = 0;
	const uint ndom = m_look->get_ndom();
	const uint nfeat = flat_features::get_nfeat();
	flat_aligner fa;
	const float open = -flat_params::m_open;
	const float ext = -flat_params::m_ext;
	fa.alloc();
	uint CurrentDomIdxT = UINT_MAX;
	for (;;)
		{
		if (m_max_secs_exceeded)
			return;
		uint dopeidx = m_NextDopeIdx++;
		if (dopeidx >= m_dope_nhit)
			{
			fa.freemem();
			return;
			}
#if SHOW_PROGRESS
		if (ThreadIdx == 0 && dopeidx%1000 == 0)
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
			doT(fa, DomIdxT);
			CurrentDomIdxT = DomIdxT;
			}
		doQ(fa, DomIdxQ, DomIdxT);
		++m_aligned_pair_count;
		}
	fa.freemem();
	}

void flat_bench::Search(
	uint ThreadCount, bool PinThreads, bool UseDope, uint MaxSecs)
	{
	Alloc();

#if SHOW_PROGRESS
	if (UseDope)
		ProgressStep(0, m_dope_nhit, "Aligning");
	else
		{
		const uint NQ = SIZE(m_Labels);
		const uint PairCount = triangle_get_K(NQ);
		ProgressStep(0, PairCount, "Aligning");
		}
#endif

	Launch(ThreadCount, PinThreads, UseDope, MaxSecs);

#if SHOW_PROGRESS
	if (UseDope)
		{
		ProgressStep(m_dope_nhit-1, m_dope_nhit, "Aligning");
		uint hits = m_ncachehits;
		uint misses = m_ncachemisses;
		ProgressLog("Cache misses %u, hits %u (%.1f%%)\n",
			misses, hits, GetPct(hits, hits+misses));
		}
	else
		{
		const uint NQ = SIZE(m_Labels);
		const uint PairCount = triangle_get_K(NQ);
		ProgressStep(PairCount-1, PairCount, "Aligning");
		}
#endif
	if (MaxSecs != UINT_MAX)
		{
		uint n = m_aligned_pair_count;
		ProgressLog("%u threads %u (%s) alignments\n",
			ThreadCount, n, IntToStr(n));
		}
	}

void flat_bench::StaticThreadBody(flat_bench *SB,
	uint ThreadIdx, bool UseDope)
	{
	if (UseDope)
		SB->ThreadBody_Dope(ThreadIdx);
	else
		SB->ThreadBody_All(ThreadIdx);
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
	vector<float> &ScalarValues)
	{
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
			|| Name == "lddtx" \
			|| Name == "lddtpow" \
			|| Name == "entropy" \
			|| Name == "rotfreetm" \
			|| Name == "revw" \
			|| StartsWith(Name, "oldts_"))
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
	flat_features::apply_weights(NameToWeight);
	}

void flat_bench::LogParams(bool show_progress) const
	{
	typedef void (*t_fn)(const char *Format, ...);
	t_fn fn = (show_progress ? ProgressLog : Log);
	fn("open=%.3g;", flat_params::m_open);
	fn("ext=%.3g;", flat_params::m_ext);
	uint nfeat = flat_features::get_nfeat();
	for (uint fi = 0; fi < nfeat; ++fi)
		fn("%s=%.3g;",
			flat_features::m_feature_names[fi].c_str(),
			flat_features::m_weights[fi]);
	fn("\n");
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
	bool need_distmxs = false;
	flat_bench::ClassifyParams(
		Names, Values, AlphaNames,
		Weights, ScalarNames, ScalarValues);

	SetScalarParams(ScalarNames, ScalarValues);

	uint n = SIZE(AlphaNames);
	asserta(SIZE(Weights) == n);
	unordered_map<string, float> NameToWeight;
	for (uint i = 0; i < n; ++i)
		{
		const string &name = AlphaNames[i];
		if (NameToWeight.find(name) != NameToWeight.end())
			Die("Dupe name in spec '%s'", name.c_str());
		NameToWeight[name] = Weights[i];
		}
	ApplyWeightsToLogOdds(NameToWeight);

	if (flat_params::need_self())
		set_selfrev_scores();
	}

void flat_bench::set_selfrev_scores()
	{
	uint ndom = m_look->get_ndom();
	if (m_self_rev_scores == 0)
		m_self_rev_scores = myalloc(float, ndom);
	flat_aligner fa;
	fa.alloc();
	for (uint domidx = 0; domidx < ndom; ++domidx)
		{
		const string &label = m_look->get_dom(domidx);
		uint L = m_fp.get_length(domidx);
		const uint8_t *prof = m_fp.get_profile(domidx);
		sid_t *distmx = 0;
		if (!m_distmxs.empty())
			distmx = m_distmxs[domidx];
		m_self_rev_scores[domidx] =
			fa.get_self_rev_score(label, prof, L);
		}
	fa.freemem();
	}

void flat_bench::set_distmxs(const vector<flat_chain_t *> &chains)
	{
	const uint M = flat_params::m_distmx_bandwidth;
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
		chaq::fill_distmx(chain, distmx);
		m_distmxs[domidx] = distmx;
		}
	}

void cmd_flat_bench()
	{
	asserta(optset_lookup);
	//asserta(optset_fapattern);
	//asserta(optset_mxpattern);
	asserta(optset_input);

	asserta(!optset_spec);
	asserta(!optset_varstr);

	const string &VarStr = g_Arg1;

	if (optset_output2) s_ftsv = CreateStdioFile(opt(output2));

// structure features
	if (optset_output3) s_f3 = CreateStdioFile(opt(output3));

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

	vector<flat_chain_t *> chains;
	read_flat_chains(opt(input), chains);
	if (optset_fapattern || optset_mxpattern)
		{
		asserta(optset_fapattern && optset_mxpattern);
		flat_features::load_alphas(feature_names, opt(mxpattern));
		FB.load_profiles(opt(fapattern));
		}
	else
		{
		flat_features::set_alphas(feature_names);
		FB.m_fp.from_chains(chains, true);//TODO hard-coded make reverse here
		}
	FB.set_distmxs(chains);

	FB.UpdateParamsFromVarStr(VarStr);
	FB.LogParams();
	FB.Alloc();
	if (flat_params::need_self())
		FB.set_selfrev_scores();
	FB.SetScalarParams(scalar_names, scalar_values);

	if (optset_label1)
		{
		asserta(optset_label2);
		flat_aligner fa;
		FB.align_pair(fa, opt(label1), opt(label2));
		return;
		}

	if (optset_thread_counts)
		{
		asserta(!optset_threads);
		const uint WARMUP_SECS = 60;
		time_t t0 = time(0);
		uint nt = GetCPUCoreCount();
		for (;;)
			{
			time_t now = time(0);
			uint elapsed_so_far = uint(now - t0);
			if (elapsed_so_far >= WARMUP_SECS)
				break;
			const uint max_secs = WARMUP_SECS - elapsed_so_far;
			ProgressLog("\nWARMUP +%u secs\n", max_secs);
			FB.Search(nt, false, optset_dope, max_secs);
			}
		ProgressLog("\nWARMUP DONE\n");
		const uint max_secs = optset_maxsecs ? opt(maxsecs) : 5;
		const uint iters = optset_iters ? opt(iters) : 5;
		vector<string> flds;
		Split(opt(thread_counts), flds, ',');
		thread_affinity ta;
		for (size_t i = 0; i < flds.size(); ++i)
			{
			const uint ThreadCount = StrToUint(flds[i]);
			bool pin = opt(no_thread_pin) ? false : ta.shouldPin(ThreadCount);
			vector<uint> alns;
			for (uint iter = 0; iter < iters; ++iter)
				{
				FB.Search(ThreadCount, pin, optset_dope, max_secs);
				alns.push_back(FB.m_aligned_pair_count);
				}
			vector<uint> order(iters);
			QuickSortOrderDesc(alns.data(), iters, order.data());
			uint median = alns[order[iters/2]];
			ProgressLog("%u threads median %u (%s) alignments pin %c\n",
				ThreadCount, median, IntToStr(median), yon(pin));
			}
		return;
		}

	uint ThreadCount = GetRequestedThreadCount();
	thread_affinity ta;
	bool pin = opt(no_thread_pin) ? false : ta.shouldPin(ThreadCount);
	FB.Search(ThreadCount, pin, optset_dope, UINT_MAX);
	FB.SetScoreOrder();
	FB.Bench();
	FB.WriteHits(opt(output), opt(include_self), opt(triangle));
	CloseStdioFile(s_ftsv);
	}
